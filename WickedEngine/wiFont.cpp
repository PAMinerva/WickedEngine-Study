#include "wiFont.h"
#include "wiRenderer.h"
#include "wiResourceManager.h"
#include "wiHelper.h"
#include "shaders/ShaderInterop_Font.h"
#include "wiBacklog.h"
#include "wiTextureHelper.h"
#include "wiRectPacker.h"
#include "wiSpinLock.h"
#include "wiPlatform.h"
#include "wiEventHandler.h"
#include "wiTimer.h"
#include "wiUnorderedMap.h"
#include "wiUnorderedSet.h"
#include "wiVector.h"
#include "wiMath.h"

#include "Utility/liberation_sans.h"
#include "Utility/stb_truetype.h"

#include <fstream>
#include <mutex>

using namespace wi::enums;
using namespace wi::graphics;

namespace wi::font
{
	namespace font_internal
	{
		enum DEPTH_TEST_MODE
		{
			DEPTH_TEST_OFF,
			DEPTH_TEST_ON,
			DEPTH_TEST_MODE_COUNT
		};
		static BlendState blendState;
		static RasterizerState rasterizerState;
		static DepthStencilState depthStencilStates[DEPTH_TEST_MODE_COUNT];

		static Shader vertexShader;
		static Shader pixelShader;
		static PipelineState PSO[DEPTH_TEST_MODE_COUNT];

		static thread_local wi::Canvas canvas;

		static Texture texture;

		struct FontStyle
		{
			std::string name;
			wi::vector<uint8_t> fontBuffer; // only used if loaded from file, need to keep alive
			stbtt_fontinfo fontInfo;
			int ascent, descent, lineGap;
			void Create(const std::string& newName, const uint8_t* data, size_t size)
			{
				name = newName;
				int offset = stbtt_GetFontOffsetForIndex(data, 0);

				if (!stbtt_InitFont(&fontInfo, data, offset))
				{
					wi::backlog::post("Failed to load font: " + name + " (file was unrecognized, it must be a .ttf file)");
				}

				stbtt_GetFontVMetrics(&fontInfo, &ascent, &descent, &lineGap);
			}
			void Create(const std::string& newName)
			{
				if (wi::helper::FileRead(newName, fontBuffer))
				{
					Create(newName, fontBuffer.data(), fontBuffer.size());
				}
				else
				{
					wi::backlog::post("Failed to load font: " + name + " (file could not be opened)");
				}
			}
		};
		static wi::vector<std::unique_ptr<FontStyle>> fontStyles;

		struct Glyph
		{
			float x;
			float y;
			float width;
			float height;
			float tc_left;
			float tc_right;
			float tc_top;
			float tc_bottom;
			const FontStyle* fontStyle = nullptr;
		};
		static wi::unordered_map<int32_t, Glyph> glyph_lookup;
		static wi::unordered_map<int32_t, wi::rectpacker::Rect> rect_lookup;
		struct Bitmap
		{
			int width;
			int height;
			int xoff;
			int yoff;
			wi::vector<uint8_t> data;
		};
		static wi::unordered_map<int32_t, Bitmap> bitmap_lookup;
		union GlyphHash
		{
			struct
			{
				uint32_t code : 16;		// character code range supported: 0 - 65535
				uint32_t height : 10;	// height supported: 0 - 1023
				uint32_t style : 5;		// number of font styles supported: 0 - 31
				uint32_t sdf : 1;		// true or false
			} bits;
			uint32_t raw;
		};
		static_assert(sizeof(GlyphHash) == sizeof(uint32_t));
		static wi::unordered_set<uint32_t> pendingGlyphs;
		static std::mutex locker;

		struct ParseStatus
		{
			Cursor cursor;
			uint32_t quadCount = 0;
			size_t last_word_begin = 0;
			bool start_new_word = false;
		};

		static thread_local wi::vector<FontVertex> vertexList;
		ParseStatus ParseText(const wchar_t* text, size_t text_length, const Params& params)
		{
			ParseStatus status;
			status.cursor = params.cursor;

			vertexList.clear();

			const float whitespace_size = (float(params.size) + params.spacingX) * 0.25f;
			const float tab_size = whitespace_size * 4;
			const float linebreak_size = (float(params.size) + params.spacingY);

			auto word_wrap = [&] {
				status.start_new_word = true;
				if (status.last_word_begin > 0 && params.h_wrap >= 0 && status.cursor.position.x >= params.h_wrap - 1)
				{
					// Word ended and wrap detected, push down last word by one line:
					const float word_offset = vertexList[status.last_word_begin].pos.x;
					for (size_t i = status.last_word_begin; i < status.quadCount * 4; ++i)
					{
						vertexList[i].pos.x -= word_offset;
						vertexList[i].pos.y += linebreak_size;
					}
					status.cursor.position.x -= word_offset;
					status.cursor.position.y += linebreak_size;
					status.cursor.size.x = std::max(status.cursor.size.x, status.cursor.position.x);
					status.cursor.size.y = std::max(status.cursor.size.y, status.cursor.position.y + linebreak_size);
				}
			};

			// update the height of the cursor size to include the next line of text
			status.cursor.size.y = status.cursor.position.y + linebreak_size;
			for (size_t i = 0; i < text_length; ++i)
			{
				int code = (int)text[i];
				GlyphHash hash;
				hash.bits.code = text[i];
				hash.bits.height = params.size;
				hash.bits.style = (uint32_t)params.style;
				hash.bits.sdf = params.isSDFRenderingEnabled() ? 1 : 0;

				if (glyph_lookup.count(hash.raw) == 0)
				{
					// glyph not packed yet, so add to pending list:
					std::scoped_lock lck(locker);
					pendingGlyphs.insert(hash.raw);
					continue;
				}

				if (code == '\n')
				{
					word_wrap();
					status.cursor.position.x = 0;
					status.cursor.position.y += linebreak_size; // \n can be in the middle of a line, so we need to add the linebreak size here
				}
				else if (code == ' ')
				{
					word_wrap();
					status.cursor.position.x += whitespace_size;
				}
				else if (code == '\t')
				{
					word_wrap();
					status.cursor.position.x += tab_size;
				}
				else if (code == '\r')
				{
				}
				else
				{
					// retrieve the user calculated glyph metrics from the lookup table
					const Glyph& glyph = glyph_lookup.at(hash.raw);
					const float glyphWidth = glyph.width;
					const float glyphHeight = glyph.height;
					const float glyphOffsetX = glyph.x;
					const float glyphOffsetY = glyph.y;
					const float fontScale = stbtt_ScaleForPixelHeight(&glyph.fontStyle->fontInfo, (float)params.size);

					const size_t vertexID = size_t(status.quadCount) * 4;
					vertexList.resize(vertexID + 4); // append additional 4 vertices to the vertex list for the current character in the text line
					status.quadCount++; // add a new quad for the current character in the text line

					if (status.start_new_word)
					{
						status.last_word_begin = vertexID;
					}
					status.start_new_word = false;

					// use the cursor (position of the prev character) to position the quad of the current character
					const float left = status.cursor.position.x + glyphOffsetX;
					const float right = left + glyphWidth;
					const float top = status.cursor.position.y + glyphOffsetY;
					const float bottom = top + glyphHeight;

					vertexList[vertexID + 0].pos = float2(left, top);
					vertexList[vertexID + 1].pos = float2(right, top);
					vertexList[vertexID + 2].pos = float2(left, bottom);
					vertexList[vertexID + 3].pos = float2(right, bottom);

					float tc_left = glyph.tc_left;
					float tc_right = glyph.tc_right;
					float tc_top = glyph.tc_top;
					float tc_bottom = glyph.tc_bottom;
					if (params.isFlippedHorizontally())
					{
						std::swap(tc_left, tc_right);
					}
					if (params.isFlippedVertically())
					{
						std::swap(tc_top, tc_bottom);
					}
					vertexList[vertexID + 0].uv = float2(tc_left, tc_top);
					vertexList[vertexID + 1].uv = float2(tc_right, tc_top);
					vertexList[vertexID + 2].uv = float2(tc_left, tc_bottom);
					vertexList[vertexID + 3].uv = float2(tc_right, tc_bottom);

					// Get the character horizontal metrics to advance the cursor position
					int advance, lsb;
					stbtt_GetCodepointHMetrics(&glyph.fontStyle->fontInfo, code, &advance, &lsb);
					status.cursor.position.x += advance * fontScale;

					status.cursor.position.x += params.spacingX;

					// If the text length is greater than 1 and the current character is not the last one,
					// calculate the kerning between the current character and the next character and apply it to the cursor position.
					// Kerning is the process of adjusting the spacing between specific pairs of characters
					// to improve the visual appearance of the text. It ensures that the spacing between
					// characters is visually pleasing and consistent. For example, the pair "AV" often
					// requires kerning to avoid excessive space between the characters.
					if (text_length > 1 && i < text_length - 1 && text[i + 1])
					{
						int code_next = (int)text[i + 1];
						int kern = stbtt_GetCodepointKernAdvance(&glyph.fontStyle->fontInfo, code, code_next);
						status.cursor.position.x += kern * fontScale;
					}
				}

				// Update the cursor size to include the next line of text.
				// The cursor size specifies the size of the entire text from the first character.
				status.cursor.size.x = std::max(status.cursor.size.x, status.cursor.position.x);
				status.cursor.size.y = std::max(status.cursor.size.y, status.cursor.position.y + linebreak_size);
			}

			word_wrap();

			return status;
		}

		thread_local static std::string char_temp_buffer;
		thread_local static std::wstring wchar_temp_buffer;
		ParseStatus ParseText(const char* text, size_t text_length, const Params& params)
		{
			// the temp buffers are used to avoid allocations of string objects:
			char_temp_buffer = text;
			wi::helper::StringConvert(char_temp_buffer, wchar_temp_buffer);
			return ParseText(wchar_temp_buffer.c_str(), wchar_temp_buffer.length(), params);
		}

		void CommitText(void* vertexList_GPU)
		{
			std::memcpy(vertexList_GPU, vertexList.data(), sizeof(FontVertex) * vertexList.size());
		}

	}
	using namespace font_internal;

	void LoadShaders()
	{
		wi::renderer::LoadShader(ShaderStage::VS, vertexShader, "fontVS.cso");
		wi::renderer::LoadShader(ShaderStage::PS, pixelShader, "fontPS.cso");

		// Create two PSOs to render the font:
		// one with no depth test, one with depth test enabled
		for (int d = 0; d < DEPTH_TEST_MODE_COUNT; ++d)
		{
			PipelineStateDesc desc;
			desc.vs = &vertexShader;
			desc.ps = &pixelShader;
			desc.bs = &blendState;
			desc.dss = &depthStencilStates[d];
			desc.rs = &rasterizerState;
			desc.pt = PrimitiveTopology::TRIANGLESTRIP;
			wi::graphics::GetDevice()->CreatePipelineState(&desc, &PSO[d]);
		}
	}
	void Initialize()
	{
		wi::Timer timer;

		// add default font if there is none yet:
		if (fontStyles.empty())
		{
			wi::vector<uint8_t> data;
			helper::Decompress(liberation_sans_zstd, sizeof(liberation_sans_zstd), data);
			// Add Liberation Sans font (as a FontStyle) to fontStyles vector and calculate font metrics: ascent, descent, lineGap
			AddFontStyle("Liberation Sans", data.data(), data.size(), true);
		}

		RasterizerState rs;
		rs.fill_mode = FillMode::SOLID;
		rs.cull_mode = CullMode::NONE;
		rs.front_counter_clockwise = true;
		rs.depth_bias = 0;
		rs.depth_bias_clamp = 0;
		rs.slope_scaled_depth_bias = 0;
		rs.depth_clip_enable = false;
		rs.multisample_enable = false;
		rs.antialiased_line_enable = false;
		rasterizerState = rs;

		BlendState bd;
		bd.render_target[0].blend_enable = true;
		bd.render_target[0].src_blend = Blend::ONE; // premultiplied blending
		bd.render_target[0].dest_blend = Blend::INV_SRC_ALPHA;
		bd.render_target[0].blend_op = BlendOp::ADD;
		bd.render_target[0].src_blend_alpha = Blend::ONE;
		bd.render_target[0].dest_blend_alpha = Blend::INV_SRC_ALPHA;
		bd.render_target[0].blend_op_alpha = BlendOp::ADD;
		bd.render_target[0].render_target_write_mask = ColorWrite::ENABLE_ALL;
		bd.independent_blend_enable = false;
		blendState = bd;

		DepthStencilState dsd;
		dsd.depth_enable = false;
		dsd.stencil_enable = false;
		depthStencilStates[DEPTH_TEST_OFF] = dsd;

		dsd.depth_enable = true;
		dsd.depth_write_mask = DepthWriteMask::ZERO;
		dsd.depth_func = ComparisonFunc::GREATER_EQUAL;
		depthStencilStates[DEPTH_TEST_ON] = dsd;

		static wi::eventhandler::Handle handle1 = wi::eventhandler::Subscribe(wi::eventhandler::EVENT_RELOAD_SHADERS, [](uint64_t userdata) { LoadShaders(); });
		LoadShaders();

		wilog("wi::font Initialized (%d ms)", (int)std::round(timer.elapsed()));
	}

	void InvalidateAtlas()
	{
		texture = {};
		glyph_lookup.clear();
		rect_lookup.clear();
		bitmap_lookup.clear();
	}
	void UpdateAtlas(float upscaling)
	{
		std::scoped_lock lck(locker);

		// upscaling used to go from logical to physical space
		upscaling = std::max(1.5f, upscaling); // add some minimum upscaling, especially for SDF
		static float upscaling_prev = 1;
		const float upscaling_rcp = 1.0f / upscaling;

		if (upscaling_prev != upscaling)
		{
			// If upscaling changed (DPI change), clear glyph caches, they will need to be re-rendered:
			InvalidateAtlas();
			upscaling_prev = upscaling;
		}

		// If there are pending glyphs, render them and repack the atlas:
		// See ParseText above
		if (!pendingGlyphs.empty())
		{
			for (int32_t raw : pendingGlyphs)
			{
				GlyphHash hash;
				hash.raw = raw;
				const int code = (int)hash.bits.code;
				const float height = (float)hash.bits.height;
				const bool is_sdf = hash.bits.sdf ? true : false;
				uint32_t style = hash.bits.style;
				FontStyle* fontStyle = fontStyles[style].get();
				int glyphIndex = stbtt_FindGlyphIndex(&fontStyle->fontInfo, code);
				if (glyphIndex == 0)
				{
					// Try fallback to an other font style that has this character:
					style = 0;
					while (glyphIndex == 0 && style < fontStyles.size())
					{
						fontStyle = fontStyles[style].get();
						glyphIndex = stbtt_FindGlyphIndex(&fontStyle->fontInfo, code);
						style++;
					}
				}

				// See stb_truetype.h
				// fontScaling can be used to scale the glyph metrics from logical to physical space, based on the font size.
				float fontScaling = stbtt_ScaleForPixelHeight(&fontStyle->fontInfo, height * upscaling);

				// xoff/yoff are the offset in pixels from the glyph origin to the top-left of the bitmap
				// yoff is the vertical offset from the glyph's origin to the top edge of the bitmap (in pixels misured with respect to the bitmap space).
				// xoff is the horizontal offset from the glyph's origin to the left edge of the bitmap (in pixels misured with respect to the bitmap space).
				// The bitmap's top-left corner is at (0, 0) in bitmap space. The glyph might not start at (0, 0) because of spacing within the character.
				// glyph's origin usually is the intersection between the baseline and the leftmost point of the glyph.
				// The space of the glyph is a logical space, with the x-axis going from left to right (similar to bitmap space) and the y-axis going
				// from bottom to top (opposite to bitmap space).
				// This means that both xoff and yoff usually are negative values because they are offsets (vectors) in bitmap space from the origin to
				// the top-left corner of the bitmap.
				Bitmap& bitmap = bitmap_lookup[hash.raw];
				bitmap.width = 0;
				bitmap.height = 0;
				bitmap.xoff = 0;
				bitmap.yoff = 0;

				if (is_sdf)
				{	// See stb_truetype.h
					unsigned char* data = stbtt_GetGlyphSDF(
						&fontStyle->fontInfo,
						fontScaling,
						glyphIndex,
						(int)SDF::padding,
						(unsigned char)SDF::onedge_value,
						SDF::pixel_dist_scale,
						&bitmap.width,
						&bitmap.height,
						&bitmap.xoff,
						&bitmap.yoff
					);
					bitmap.data.resize(bitmap.width * bitmap.height);
					if (data) std::memcpy(bitmap.data.data(), data, bitmap.data.size());
					stbtt_FreeSDF(data, nullptr);
				}
				else
				{
					unsigned char* data = stbtt_GetGlyphBitmap(
						&fontStyle->fontInfo,
						fontScaling,
						fontScaling,
						glyphIndex,
						&bitmap.width,
						&bitmap.height,
						&bitmap.xoff,
						&bitmap.yoff
					);
					bitmap.data.resize(bitmap.width * bitmap.height);
					if (data) std::memcpy(bitmap.data.data(), data, bitmap.data.size());
					stbtt_FreeBitmap(data, nullptr);
				}

				wi::rectpacker::Rect rect = {};
				rect.w = bitmap.width + 2;
				rect.h = bitmap.height + 2;
				rect.id = hash.raw;
				rect_lookup[hash.raw] = rect;

				Glyph& glyph = glyph_lookup[hash.raw];
				glyph.x = float(bitmap.xoff) * upscaling_rcp; // see PhysicalToLogical in wiCanvas.h
				// ascent * fontScaling is the offset from the the top of the glyph to the baseline,
				// in pixels and scaled based on the font size.
				// By adding yoff to ascent * fontScaling, we get the offset from the top of the bitmap
				// to the top of the glyph, misured in pixels with respect to the bitmap space (similar to xoff,
				// which is the offset from the left of the glyph to the left of the bitmap).
				// glyph.x is the logical offset from the leftmost point of the glyph to the left border of the glyph
				// in the glyph space.
				// glyph.y is the logical offset from the top border of the glyph to the uppermost point of the glyph
				// in the glyph space.
				glyph.y = (float(bitmap.yoff) + float(fontStyle->ascent) * fontScaling) * upscaling_rcp;
				glyph.width = float(bitmap.width) * upscaling_rcp;   // see PhysicalToLogical in wiCanvas.h
				glyph.height = float(bitmap.height) * upscaling_rcp; // see PhysicalToLogical in wiCanvas.h
				glyph.fontStyle = fontStyle;
			}
			pendingGlyphs.clear();

			// Setup packer, this will allocate memory if needed:
			static thread_local wi::rectpacker::State packer;
			packer.clear();
			for (auto& it : rect_lookup)
			{
				packer.add_rect(it.second);
			}

			// Perform packing and process the result if successful:
			// pack multiple characters into one atlas
			if (packer.pack(4096))
			{
				// Retrieve texture atlas dimensions:
				const int atlasWidth = packer.width;
				const int atlasHeight = packer.height;
				const float inv_width = 1.0f / atlasWidth;
				const float inv_height = 1.0f / atlasHeight;

				// Create the CPU-side texture atlas and fill with transparency (0):
				// so each texel is an 8-bit value where 0 is no coverage (transparent), 255 is fully covered (opaque).
				wi::vector<uint8_t> atlas(size_t(atlasWidth) * size_t(atlasHeight));
				std::fill(atlas.begin(), atlas.end(), 0);

				// Iterate all packed glyph rectangles:
				for (auto& rect : packer.rects)
				{
					// Above the rect is built like this:
					// 
					// wi::rectpacker::Rect rect = {};
					// rect.w = bitmap.width + 2;
					// rect.h = bitmap.height + 2;
					// 
					// so this means that, for example, if the bitmap were 8x8, the rect will
					// be 10x10 and packer.pack will calculate x and y based on the 10x10 rect.
					// Now, we want to adjust the rect to fit the actual bitmap size, but
					// subtracting 2 from the width and height of the rect is not enough,
					// because this operation will move the 8x8 rect to the right and down
					// by 2 pixels in the 10x10 area of the atlas reserved for the glyph.
					// So, we need to move the rect to the left and up by 1 pixel to make
					// the 8x8 bitmap fit the 10x10 rect in the center.
					rect.x += 1;
					rect.y += 1;
					rect.w -= 2;
					rect.h -= 2;

					const int32_t hash = rect.id;
					//const wchar_t code = codefromhash(hash);
					//const int style = stylefromhash(hash);
					//const float height = (float)heightfromhash(hash);
					Glyph& glyph = glyph_lookup[hash];
					Bitmap& bitmap = bitmap_lookup[hash];

					// Copy, row by row, the glyph bitmap to the CPU-side texture atlas
					// Note that a row in the glyph bitmap is a row in the CPU-side texture atlas but with a different pitch
					// since the atlas includes all the rows of all the glyphs.
					for (int row = 0; row < bitmap.height; ++row)
					{
 						uint8_t* dst = atlas.data() + rect.x + (rect.y + row) * atlasWidth;
						uint8_t* src = bitmap.data.data() + row * bitmap.width;
						std::memcpy(dst, src, bitmap.width);
					}

					// Compute texture coordinates for the glyph:
					glyph.tc_left = float(rect.x);
					glyph.tc_right = glyph.tc_left + float(rect.w);
					glyph.tc_top = float(rect.y);
					glyph.tc_bottom = glyph.tc_top + float(rect.h);

					// Normalize texture coordinates:
					glyph.tc_left *= inv_width;
					glyph.tc_right *= inv_width;
					glyph.tc_top *= inv_height;
					glyph.tc_bottom *= inv_height;
				}

				// Upload the CPU-side texture atlas bitmap to the GPU:
				// The texture atlas is a 2D texture where each texel is 8-bit (R8_UNORM) with 0 is no coverage (transparent),
				// 255 is fully covered (opaque).
				// Also create an SRV to the texture atlas and write it to an appropriate heap (copied in a bindless one as well if
				// slots are available; if that's the case, an index to the descriptor is stored in the internal state of the texture).
				wi::texturehelper::CreateTexture(texture, atlas.data(), atlasWidth, atlasHeight, Format::R8_UNORM);
				GetDevice()->SetName(&texture, "wi::font::texture");
			}
			else
			{
				assert(0); // rect packing failure
			}
		}

	}
	const Texture* GetAtlas()
	{
		return &texture;
	}
	int AddFontStyle(const std::string& fontName)
	{
		std::scoped_lock lck(locker);
		for (size_t i = 0; i < fontStyles.size(); i++)
		{
			const FontStyle& fontStyle = *fontStyles[i];
			if (fontStyle.name.compare(fontName) == 0)
			{
				return int(i);
			}
		}
		fontStyles.push_back(std::make_unique<FontStyle>());
		fontStyles.back()->Create(fontName);
		InvalidateAtlas(); // invalidate atlas, in case there were missing glyphs, upon adding new font style they could become valid
		return int(fontStyles.size() - 1);
	}
	int AddFontStyle(const std::string& fontName, const uint8_t* data, size_t size, bool copyData)
	{
		std::scoped_lock lck(locker);
		for (size_t i = 0; i < fontStyles.size(); i++)
		{
			const FontStyle& fontStyle = *fontStyles[i];
			if (fontStyle.name.compare(fontName) == 0)
			{
				return int(i);
			}
		}
		fontStyles.push_back(std::make_unique<FontStyle>());
		if (copyData)
		{
			fontStyles.back()->fontBuffer.resize(size);
			std::memcpy(fontStyles.back()->fontBuffer.data(), data, size);
			data = fontStyles.back()->fontBuffer.data();
		}
		fontStyles.back()->Create(fontName, data, size);
		InvalidateAtlas(); // invalidate atlas, in case there were missing glyphs, upon adding new font style they could become valid
		return int(fontStyles.size() - 1);
	}

	template<typename T>
	Cursor Draw_internal(const T* text, size_t text_length, const Params& params, CommandList cmd)
	{
		if (text_length <= 0)
		{
			return Cursor();
		}
		// Parse the text line to get the number of quads and the cursor position and size.
		// The cursor is returned to be used as the starting point for the next text line.
		// ParseText will also build the vertex buffer for the text line by using the current cursor position.
		ParseStatus status = ParseText(text, text_length, params);

		// quadCount is the number of quads, each wrapped by a texture representing a text character
		if (status.quadCount > 0)
		{
			GraphicsDevice* device = wi::graphics::GetDevice();
			// Allocate enough upload memory to store the vertex buffer with all the vertices composing the quads for a text line.
			// All text lines will be stored in the same buffer, so the allocation needs to be built iteratively (see AllocateGPU).
			// A reference to the shared buffer containing all text lines can be retrieved from the internal state of the
			// command list associated with the current frame.
			// Also create an SRV to shared buffer and write it to an appropriate heap (copied in a bindless one as well if
			// slots are available; if that's the case, an index to the descriptor is stored in the internal state of the texture).
			// Remember that the index of the descriptor (in the bindless heap or in the subresources_srv array) can be
			// retrieved from the internal state of the buffer.
			GraphicsDevice::GPUAllocation mem = device->AllocateGPU(sizeof(FontVertex) * status.quadCount * 4, cmd);
			if (!mem.IsValid())
			{
				// Returning the current cursor, the next time the function is called, we can
				// continue rendering the text from where it was interrupted.
				return status.cursor;
			}
			CommitText(mem.data); // copy the vertex data to the upload memory

			FontConstants font = {};
			font.buffer_index = device->GetDescriptorIndex(&mem.buffer, SubresourceType::SRV); // get the index of the SRV (describing the vertex buffer) in the heap
			font.buffer_offset = (uint32_t)mem.offset; // offset of the vertex buffer portion containing the current text line
			font.texture_index = device->GetDescriptorIndex(&texture, SubresourceType::SRV); // get the index of the SRV (describing the texture atlas) in the heap
			if (font.buffer_index < 0 || font.texture_index < 0) // does it support bindless only?
			{
				return status.cursor;
			}

			device->EventBegin("Font", cmd); // Starts a user-defined event for a timing capture of CPU activity, to be displayed in PIX

			// If the PSO has not been created yet, save the PipelineState passed as first param as the active one in the command list and
			// set it dirty so that a PSO based on it can be created before the actual draw call.
			// It also set the root signature and invalidates all root bindings if necessary.
			device->BindPipelineState(&PSO[params.isDepthTestEnabled()], cmd);

			using namespace wi::math;
			XMFLOAT4 color = XMFLOAT4(1, 1, 1, 1);
			float softness = 0;
			float bolden = 0;
			float hdr_scaling = 1;
			uint32_t flags = 0;

			if (params.isSDFRenderingEnabled())
			{
				flags |= FONT_FLAG_SDF_RENDERING;
			}
			if (params.isHDR10OutputMappingEnabled())
			{
				flags |= FONT_FLAG_OUTPUT_COLOR_SPACE_HDR10_ST2084;
			}
			if (params.isLinearOutputMappingEnabled())
			{
				flags |= FONT_FLAG_OUTPUT_COLOR_SPACE_LINEAR;
				hdr_scaling = params.hdr_scaling;
			}

			XMFLOAT3 offset = XMFLOAT3(0, 0, 0);
			float vertical_flip = params.customProjection == nullptr ? 1.0f : -1.0f;
			if (params.h_align == WIFALIGN_CENTER)
				offset.x -= status.cursor.size.x / 2;
			else if (params.h_align == WIFALIGN_RIGHT)
				offset.x -= status.cursor.size.x;
			if (params.v_align == WIFALIGN_CENTER)
				offset.y -= status.cursor.size.y / 2 * vertical_flip;
			else if (params.v_align == WIFALIGN_BOTTOM)
				offset.y -= status.cursor.size.y * vertical_flip;

			XMMATRIX M = XMMatrixTranslation(offset.x, offset.y, offset.z);
			M = M * XMMatrixScaling(params.scaling, params.scaling, params.scaling);
			M = M * XMMatrixRotationZ(params.rotation);

			if (params.customRotation != nullptr)
			{
				M = M * (*params.customRotation);
			}

			M = M * XMMatrixTranslation(params.position.x, params.position.y, params.position.z);

			if (params.customProjection != nullptr)
			{
				M = XMMatrixScaling(1, -1, 1) * M; // reason: screen projection is Y down (like UV-space) and that is the common case for image rendering. But custom projections will use the "world space"
				M = M * (*params.customProjection);
			}
			else
			{
				// Asserts will check that a proper canvas was set for this cmd with wi::image::SetCanvas()
				//	The canvas must be set to have dpi aware rendering
				assert(canvas.width > 0);
				assert(canvas.height > 0);
				assert(canvas.dpi > 0);
				M = M * canvas.GetProjection();
			}

			if (params.shadowColor.getA() > 0)
			{
				// font shadow render:
				XMStoreFloat4x4(&font.transform, XMMatrixTranslation(params.shadow_offset_x, params.shadow_offset_y, 0) * M);
				color = params.shadowColor;
				color.x *= params.shadow_intensity;
				color.y *= params.shadow_intensity;
				color.z *= params.shadow_intensity;
				font.color = pack_half4(color);
				bolden = params.shadow_bolden;
				softness = params.shadow_softness * 0.5f;
				font.softness_bolden_hdrscaling = pack_half3(softness, bolden, hdr_scaling);
				font.softness_bolden_hdrscaling.y |= flags << 16u;
				device->BindDynamicConstantBuffer(font, CBSLOT_FONT, cmd);

				device->DrawInstanced(4, status.quadCount, 0, 0, cmd);
			}

			// font base render:
			XMStoreFloat4x4(&font.transform, M);
			color = params.color;
			color.x *= params.intensity;
			color.y *= params.intensity;
			color.z *= params.intensity;
			font.color = pack_half4(color);
			bolden = params.bolden;
			softness = params.softness * 0.5f;
			font.softness_bolden_hdrscaling = pack_half3(softness, bolden, hdr_scaling);
			font.softness_bolden_hdrscaling.y |= flags << 16u;
			// The font constants for all text lines will be also loaded in the same buffer containing the vertex buffer with all text lines.
			// However, the portion of buffer containing the constants for the current text line will be accessed in the
			// vertex shader through a constant buffer (slot b0).
			// The binder of the command list will be updated to store a reference to the buffer containing the font constants for all text lines
			// and the offset to the font constants for the current text line.
			// In the binder the related root parameter will be marked as dirty so that we can update the related root argument before the draw call.
			device->BindDynamicConstantBuffer(font, CBSLOT_FONT, cmd);

			// Check if the PSO needs to be created and set before invoking the actual draw call.
			// Use instanced rendering to draw the four vertices of a quad quadCount times, so that the entire text line is rendered.
			// The vertex data will be retrieved in the VS from the bindless_buffers array
			// and the FontConstants data specified above from the ConstantBuffer font, see fontVS.hlsl, globals.hlsli, ShaderInterop_Font.h and ShaderInterop.h 
			device->DrawInstanced(4, status.quadCount, 0, 0, cmd);

			device->EventEnd(cmd);
		}

		return status.cursor;
	}

	void SetCanvas(const wi::Canvas& current_canvas)
	{
		canvas = current_canvas;
	}

	Cursor Draw(const char* text, size_t text_length, const Params& params, CommandList cmd)
	{
		return Draw_internal(text, text_length, params, cmd);
	}
	Cursor Draw(const wchar_t* text, size_t text_length, const Params& params, CommandList cmd)
	{
		return Draw_internal(text, text_length, params, cmd);
	}
	Cursor Draw(const char* text, const Params& params, CommandList cmd)
	{
		return Draw_internal(text, strlen(text), params, cmd);
	}
	Cursor Draw(const wchar_t* text, const Params& params, CommandList cmd)
	{
		return Draw_internal(text, wcslen(text), params, cmd);
	}
	Cursor Draw(const std::string& text, const Params& params, CommandList cmd)
	{
		return Draw_internal(text.c_str(), text.length(), params, cmd);
	}
	Cursor Draw(const std::wstring& text, const Params& params, CommandList cmd)
	{
		return Draw_internal(text.c_str(), text.length(), params, cmd);
	}

	XMFLOAT2 TextSize(const char* text, size_t text_length, const Params& params)
	{
		if (text_length == 0)
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text, text_length, params).cursor.size;
	}
	XMFLOAT2 TextSize(const wchar_t* text, size_t text_length, const Params& params)
	{
		if (text_length == 0)
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text, text_length, params).cursor.size;
	}
	XMFLOAT2 TextSize(const char* text, const Params& params)
	{
		size_t text_length = strlen(text);
		if (text_length == 0)
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text, text_length, params).cursor.size;
	}
	XMFLOAT2 TextSize(const wchar_t* text, const Params& params)
	{
		size_t text_length = wcslen(text);
		if (text_length == 0)
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text, text_length, params).cursor.size;
	}
	XMFLOAT2 TextSize(const std::string& text, const Params& params)
	{
		if (text.empty())
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text.c_str(), text.length(), params).cursor.size;
	}
	XMFLOAT2 TextSize(const std::wstring& text, const Params& params)
	{
		if (text.empty())
		{
			return XMFLOAT2(0, 0);
		}
		return ParseText(text.c_str(), text.length(), params).cursor.size;
	}

	Cursor TextCursor(const char* text, size_t text_length, const Params& params)
	{
		if (text_length == 0)
		{
			return {};
		}
		return ParseText(text, text_length, params).cursor;
	}
	Cursor TextCursor(const wchar_t* text, size_t text_length, const Params& params)
	{
		if (text_length == 0)
		{
			return {};
		}
		return ParseText(text, text_length, params).cursor;
	}
	Cursor TextCursor(const char* text, const Params& params)
	{
		size_t text_length = strlen(text);
		if (text_length == 0)
		{
			return {};
		}
		return ParseText(text, text_length, params).cursor;
	}
	Cursor TextCursor(const wchar_t* text, const Params& params)
	{
		size_t text_length = wcslen(text);
		if (text_length == 0)
		{
			return {};
		}
		return ParseText(text, text_length, params).cursor;
	}
	Cursor TextCursor(const std::string& text, const Params& params)
	{
		if (text.empty())
		{
			return {};
		}
		return ParseText(text.c_str(), text.length(), params).cursor;
	}
	Cursor TextCursor(const std::wstring& text, const Params& params)
	{
		if (text.empty())
		{
			return {};
		}
		return ParseText(text.c_str(), text.length(), params).cursor;
	}

	float TextWidth(const char* text, size_t text_length, const Params& params)
	{
		return TextSize(text, text_length, params).x;
	}
	float TextWidth(const wchar_t* text, size_t text_length, const Params& params)
	{
		return TextSize(text, text_length, params).x;
	}
	float TextWidth(const char* text, const Params& params)
	{
		return TextSize(text, params).x;
	}
	float TextWidth(const wchar_t* text, const Params& params)
	{
		return TextSize(text, params).x;
	}
	float TextWidth(const std::string& text, const Params& params)
	{
		return TextSize(text, params).x;
	}
	float TextWidth(const std::wstring& text, const Params& params)
	{
		return TextSize(text, params).x;
	}

	float TextHeight(const char* text, size_t text_length, const Params& params)
	{
		return TextSize(text, text_length, params).y;
	}
	float TextHeight(const wchar_t* text, size_t text_length, const Params& params)
	{
		return TextSize(text, text_length, params).y;
	}
	float TextHeight(const char* text, const Params& params)
	{
		return TextSize(text, params).y;
	}
	float TextHeight(const wchar_t* text, const Params& params)
	{
		return TextSize(text, params).y;
	}
	float TextHeight(const std::string& text, const Params& params)
	{
		return TextSize(text, params).y;
	}
	float TextHeight(const std::wstring& text, const Params& params)
	{
		return TextSize(text, params).y;
	}

}
