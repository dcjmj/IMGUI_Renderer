//Render Buffer
#include "rendering_API.h"
#include "lodepng.h"
#include "cyPoint.h"

void saveTexture2DAsImage(GLuint texture, const std::string& filename) {
	std::cout << "now is " << texture << std::endl;
	glBindTexture(GL_TEXTURE_2D, texture);

	// Get the width and height of the texture
	GLint width, height;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);

	// Read back the texture data
	std::vector<float> outputData(width * height * 4);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, outputData.data());

	std::vector<unsigned char> pixelData(width * height * 4);

	for (size_t i = 0; i < outputData.size(); ++i) {
		pixelData[i] = static_cast<unsigned char>(std::clamp(outputData[i] * 255.0f, 0.0f, 255.0f));
	}

	// Save as PNG
	unsigned error = lodepng::encode(filename, pixelData, width, height);
	if (error) {
		std::cerr << "PNG encoder error: " << lodepng_error_text(error) << std::endl;
	}
}
// -------------------------------------------------------------------------- Physically-based Rendering ---------------------------------------------------------------------
// GL related setup
void OpenGLRenderBuffer::Initialize(bool useDepthBuffer, int num_colorbuffer)
{
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevBuffer);

	glGenFramebuffers(1, &buffer);
	glBindFramebuffer(GL_FRAMEBUFFER, buffer);

	glGenTextures(1, &bufferTex);
	glBindTexture(GL_TEXTURE_2D, bufferTex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	if (useDepthBuffer) {
		glGenRenderbuffers(1, &bufferDepth);
		glBindRenderbuffer(GL_RENDERBUFFER, bufferDepth);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, bufferDepth);
	}

	if (num_colorbuffer > 1) {
		glGenTextures(1, &buffer2Tex);
		glBindTexture(GL_TEXTURE_2D, buffer2Tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, buffer2Tex, 0);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, bufferTex, 0);

		GLenum DrawBuffers[2] = { GL_COLOR_ATTACHMENT0,  GL_COLOR_ATTACHMENT1 };
		glDrawBuffers(2, DrawBuffers);
	}
	else {
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, bufferTex, 0);

		GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
		glDrawBuffers(1, DrawBuffers);
		buffer2Tex = -1;
	}

	int bufferReady = glIsFramebuffer(buffer);
	assert(bufferReady);

	glBindFramebuffer(GL_FRAMEBUFFER, prevBuffer);
}

//-------------------------------------------------------------------------------

void OpenGLRenderBuffer::Resize(int width, int height, Format format, int numChannels)
{
	assert(numChannels <= 4 && numChannels > 0);

	const static GLint internalFormat[] = {
		GL_R8I,   GL_RG8I,   GL_RGB8I,   GL_RGBA8I,
		GL_R16I,  GL_RG16I,  GL_RGB16I,  GL_RGBA16I,
		GL_R32I,  GL_RG32I,  GL_RGB32I,  GL_RGBA32I,
		GL_R8UI,  GL_RG8UI,  GL_RGB8UI,  GL_RGBA8UI,
		GL_R16UI, GL_RG16UI, GL_RGB16UI, GL_RGBA16UI,
		GL_R32UI, GL_RG32UI, GL_RGB32UI, GL_RGBA32UI,
		GL_R16F,  GL_RG16F,  GL_RGB16F,  GL_RGBA16F,
		GL_R32F,  GL_RG32F,  GL_RGB32F,  GL_RGBA32F,
	};
	const static GLenum gformat[] = {
		GL_RED,         GL_RG,         GL_RGB,         GL_RGBA,
		GL_RED_INTEGER, GL_RG_INTEGER, GL_RGB_INTEGER, GL_RGBA_INTEGER,
		GL_RED_INTEGER, GL_RG_INTEGER, GL_RGB_INTEGER, GL_RGBA_INTEGER,
		GL_RED,         GL_RG,         GL_RGB,         GL_RGBA,
		GL_RED_INTEGER, GL_RG_INTEGER, GL_RGB_INTEGER, GL_RGBA_INTEGER,
		GL_RED_INTEGER, GL_RG_INTEGER, GL_RGB_INTEGER, GL_RGBA_INTEGER,
		GL_RED,         GL_RG,         GL_RGB,         GL_RGBA,
		GL_RED,         GL_RG,         GL_RGB,         GL_RGBA,
	};
	const static GLenum type[] = {
		GL_BYTE,
		GL_SHORT,
		GL_INT,
		GL_UNSIGNED_BYTE,
		GL_UNSIGNED_SHORT,
		GL_UNSIGNED_INT,
		GL_FLOAT,
		GL_FLOAT,
	};

	int fi = 4 * format + numChannels - 1;

	glBindTexture(GL_TEXTURE_2D, bufferTex);
	glTexImage2D(GL_TEXTURE_2D, 0, internalFormat[fi], width, height, 0, gformat[fi], type[format], 0);

	if (buffer2Tex >= 0) {
		glBindTexture(GL_TEXTURE_2D, buffer2Tex);
		glTexImage2D(GL_TEXTURE_2D, 0, internalFormat[fi], width, height, 0, gformat[fi], type[format], 0);
	}

	if (bufferDepth) {
		glBindRenderbuffer(GL_RENDERBUFFER, bufferDepth);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
	}

	Bind();
	bool bufferComplete = (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
	assert(bufferComplete);
	Unbind();
}

//-------------------------------------------------------------------------------

void OpenGLRenderBuffer::Delete()
{
	if (bufferDepth != 0) {
		glDeleteBuffers(1, &bufferDepth);
		bufferDepth = 0;
	}
	if (buffer != 0) {
		glDeleteTextures(1, &bufferTex);
		glDeleteBuffers(1, &buffer);
		buffer = 0;
	}
	if (buffer2Tex != -1) {
		glDeleteTextures(1, &bufferTex);
		glDeleteBuffers(1, &buffer);
		buffer = 0;
	}
}

//-------------------------------------------------------------------------------

void OpenGLRenderBuffer::Bind()
{
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, buffer);
}

//-------------------------------------------------------------------------------

void OpenGLRenderBuffer::Unbind()
{
	glBindFramebuffer(GL_FRAMEBUFFER, prevBuffer);
}

void GLWidget::UpdateConfig2OfflineYarn(Fiber::Yarn* y)
{
	y->plys.resize((int)fiberData.g_ply_num);

	y->use_flyaways = fiberData.g_use_flyaways;
	y->z_step_size = fiberData.g_z_step_size;
	y->z_step_num = fiberData.g_z_step_num;
	y->clock_wise = fiberData.g_yarn_clock_wise;
	y->use_migration = fiberData.g_use_migration;
	y->yarn_alpha = fiberData.g_yarn_alpha;
	y->yarn_radius = fiberData.g_yarn_radius;

	for (int i = 0; i < (int)y->plys.size(); i++) {
		Fiber::Ply& p = y->plys[i];
		p.fibers.resize((int)fiberData.g_fiber_num);
		p.flyaway_hair_density = fiberData.g_flyaway_hair_density;
		p.flyaway_hair_ze_mu = fiberData.g_flyaway_hair_ze_mu;
		p.flyaway_hair_ze_sigma = fiberData.g_flyaway_hair_ze_sigma;
		p.flyaway_hair_r0_mu = fiberData.g_flyaway_hair_r0_mu;
		p.flyaway_hair_r0_sigma = fiberData.g_flyaway_hair_r0_sigma;
		p.flyaway_hair_re_mu = fiberData.g_flyaway_hair_re_mu;
		p.flyaway_hair_re_sigma = fiberData.g_flyaway_hair_re_sigma;
		p.flyaway_hair_pe_mu = fiberData.g_flyaway_hair_pe_mu;
		p.flyaway_hair_pe_sigma = fiberData.g_flyaway_hair_pe_sigma;
		p.flyaway_loop_density = fiberData.g_flyaway_loop_density;
		p.flyaway_loop_r1_mu = fiberData.g_flyaway_loop_r1_mu;
		p.flyaway_loop_r1_sigma = fiberData.g_flyaway_loop_r1_sigma;
		p.clock_wise = fiberData.g_fiber_clock_wise;
		p.epsilon = fiberData.g_epsilon;
		p.R_max = fiberData.g_r_max;
		p.beta = fiberData.g_beta;
		p.alpha = fiberData.g_alpha;
		p.s_i = fiberData.g_s_i;
		p.rho_min = fiberData.g_rho_min;
		p.rho_max = fiberData.g_rho_max;
		p.ellipse_long = fiberData.g_ellipse_long;
		p.ellipse_short = fiberData.g_ellipse_short;
	}

	for (int i = 0; i < 3; i++) {
		y->aabb_micro_ct.pMin[i] = fiberData.g_aabb_micro_ct_pMin[i];
		y->aabb_micro_ct.pMax[i] = fiberData.g_aabb_micro_ct_pMax[i];
	}
}

void GLWidget::InitOfflineYarns()
{
	//////////////////////////////////////////////////////////////////////////
	// create offline ply buffer

	// generate plys

	offlinePlys = new Fiber::Yarn();
	UpdateConfig2OfflineYarn(offlinePlys);
	offlinePlys->fiber_radius = 0.158f;
	offlinePlys->is_uniform = false;
	offlinePlys->simulate_ply();

	// create cross section samples for real-time rendering
	offlinePlys->createCrossSectionSamples();

	// roll plys
	const std::string gen_fibers_file = "gen_fibers.txt";
	offlineYarns = new Fiber::Yarn();

	// use previous loaded config to offline yarn generator
	UpdateConfig2OfflineYarn(offlineYarns);
	offlineYarns->fiber_radius = 0.158f;
	offlineYarns->z_step_num = 38;
	offlineYarns->aabb_micro_ct.pMin.z() = -0.5f * offlineYarns->z_step_num * offlineYarns->z_step_size;
	offlineYarns->aabb_micro_ct.pMax.z() = (0.5f * offlineYarns->z_step_num - 1) * offlineYarns->z_step_size;
	offlineYarns->use_flyaways = false;
	offlineYarns->yarn_radius = 0.02866;
	offlineYarns->simulate_fly_away();

	// delete previous allocated buffer
	for (int i = 0; i < (int)offlineYarnVertices.size(); i++)
		glDeleteBuffers(1, &offlineYarnBuffer[i]);
	//delete[] offlineYarnBuffer;

	for (int i = 0; i < (int)offlineYarnVertices.size(); i++)
		offlineYarnVertices[i].clear();
	offlineYarnVertices.clear();

	for (int ply = 0; ply < offlineYarns->plys.size(); ply++) {
		int num_hair_fibers = offlineYarns->plys[ply].num_hair_fibers;
		int fiber_num = (int)offlineYarns->plys[ply].fibers.size();
		for (int f = 0; f < fiber_num; f++) {
			Fiber::Fiber& fiber = offlineYarns->plys[ply].fibers[f];
			std::vector<ks::vec3> fiberVertices;
			int fiber_vertex_num = (int)fiber.vertices.size();
			for (int v = 0; v < fiber_vertex_num; v++) {
				fiberVertices.push_back(ks::vec3(fiber.vertices[v].z(), fiber.vertices[v].y(), fiber.vertices[v].x()));
			}

			offlineYarnVertices.push_back(fiberVertices);
		}
	}

	// put final offline yarn data to buffer
	offlineYarnBuffer = new GLuint[(int)offlineYarnVertices.size()];
	for (int i = 0; i < (int)offlineYarnVertices.size(); i++)
	{
		glGenBuffers(1, &offlineYarnBuffer[i]);
		glBindBuffer(GL_ARRAY_BUFFER, offlineYarnBuffer[i]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(ks::vec3) * (offlineYarnVertices[i].size()), &offlineYarnVertices[i][0], GL_STATIC_DRAW);
	}
}

void GLWidget::Init1DCylinderNormTexture()
{
	std::cout << "Loading 1D cylinder norm texture... ";
	glGenTextures(1, &cylinderTex);

	int num = 100;
	std::vector<float>	cylinderTextureData;
	cylinderTextureData.resize(num);

	for (int i = 0; i < num; i++)
		//cylinderTextureData[i] = float(i) / (num);
		cylinderTextureData[i] = sqrt(1.0f - (float(i) / (num)-0.5f) * (float(i) / (num)-0.5f) * 4.0f);

	glBindTexture(GL_TEXTURE_1D, cylinderTex);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_R32F, num, 0, GL_RED, GL_FLOAT, cylinderTextureData.data());
	glGenerateMipmap(GL_TEXTURE_1D);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

//Loading my Greyscale image for GIAO 
bool loadPNG(const std::string& filename, std::vector<unsigned char>& image, unsigned& width, unsigned& height) {
	lodepng::State state;
	unsigned error = lodepng::decode(image, width, height, filename, LCT_GREY, 8);
	if (error) {
		std::cerr << "Error loading PNG: " << lodepng_error_text(error) << std::endl;
		return false;
	}
	return true;
}

//Loading my Greyscale image for GIAO 
bool loadRGBPNG(const std::string& filename, std::vector<unsigned char>& image, unsigned& width, unsigned& height) {
	lodepng::State state;
	unsigned error = lodepng::decode(image, width, height, filename, LCT_RGB, 8);
	if (error) {
		std::cerr << "Error loading PNG: " << lodepng_error_text(error) << std::endl;
		return false;
	}
	return true;
}

float gammaCorrection(float value, float gamma) {
	if (value < 0.0f) {
		return 0.0f;
	}

	return std::pow(value, 1.0f / gamma);
}

void GLWidget::InitCoreAOTexture() {
	std::cout << "Read CoreAO texture: ";

	std::vector<std::string> filenames = {
	   "./high_res_AO/AO_0.png",
	   "./high_res_AO/AO_1.png",
	   "./high_res_AO/AO_2.png",
	   "./high_res_AO/AO_3.png",
	   "./high_res_AO/AO_4.png",
	   "./high_res_AO/AO_5.png",
	   "./high_res_AO/AO_6.png",
	   "./high_res_AO/AO_7.png"
	};

	std::vector<std::vector<float>> images;
	unsigned width = 0, height = 0;
	unsigned xmin_R, ymin_R;

	int num_angle = 0;

	std::vector<unsigned char> image_z;
	unsigned img_width_z = 0, img_height_z = 0;
	std::string z_filename = "D:/StitchMesh4.0/high_res_AO/offset_z.png";
	if (!loadRGBPNG(z_filename, image_z, img_width_z, img_height_z)) {
		std::cerr << "Failed to load image: " << z_filename << std::endl;
		return;
	}

	for (const auto& filename : filenames) {
		std::vector<unsigned char> image;
		unsigned img_width = 0, img_height = 0;

		if (!loadRGBPNG(filename, image, img_width, img_height)) {
			std::cerr << "Failed to load image: " << filename << std::endl;
			return;
		}

		GLfloat* pixels = new GLfloat[img_width * img_height * 4];
		for (size_t i = 0; i < image.size(); ++i) {
			pixels[i / 3 * 4 + i % 3] = gammaCorrection(image[i] / 255.0f, 1.f / 1);
			if (i % 3 == 2) pixels[i / 3 * 4 + 3] = gammaCorrection(image_z[i - 2] / 255.0f, 1.f / 1);
		}


		//get the empty pixel out of origin image
		// find min max
		int xmin = img_width, ymin = img_height;
		int xmax = 0, ymax = 0;
		float depth_max = 0.0;
		float depth_min = 255.0;
		for (int y = 0; y < img_height; y++)
			for (int x = 0; x < img_width; x++)
			{
				int startAddressOfPixel = ((y * img_width) + x) * 4;
				if ((pixels[startAddressOfPixel + 3] > 0.001))
				{
					if (xmin > x) xmin = x;
					if (ymin > y) ymin = y;
					if (xmax < x) xmax = x;
					if (ymax < y) ymax = y;
					depth_max = (depth_max < pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_max;
					depth_min = (depth_min > pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_min;
					//printf("pixel at 
				}
			}

		std::cout << "Core AO texture bounding box: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
		std::cout << "AO ILLUMINANCE RANGE: " << depth_max << " " << depth_min << std::endl;

		// clamp the border
		ymax += 10;

		std::cout << "height: " << ymax - ymin + 1 << " width: " << xmax - xmin + 1 << std::endl;


		// create core texture
		g_core_AO_texture_width = xmax - xmin + 1;
		g_core_AO_texture_height = ymax - ymin + 1;
		float convert_depth_max = 0.0;
		float convert_depth_min = 255.0;
		std::vector<GLfloat> corePixels(g_core_AO_texture_width * g_core_AO_texture_height * 4);
		float r_min = 255.0;
		float r_max = -255.0;
		float g_min = 255.0;
		float g_max = -255.0;
		float b_min = 255.0;
		float b_max = -255.0;
		for (int y = ymin; y <= ymax - 10; y++)
		{
			for (int x = xmin; x <= xmax; x++)
			{
				int startAddressOfPixel = ((y * img_width) + x) * 4;
				int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
				corePixels[startAddressOfCorePixel + 0] = pixels[startAddressOfPixel + 0];
				corePixels[startAddressOfCorePixel + 1] = pixels[startAddressOfPixel + 1];
				corePixels[startAddressOfCorePixel + 2] = pixels[startAddressOfPixel + 2];
				if (pixels[startAddressOfPixel + 3] < 0.001)
				{
					corePixels[startAddressOfCorePixel + 3] = 0.0;		// black part
				}
				else
				{
					r_min = (r_min > pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_min;
					g_min = (g_min > pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_min;
					b_min = (b_min > pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_min;

					r_max = (r_min < pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_max;
					g_max = (g_min < pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_max;
					b_max = (b_min < pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_max;

					corePixels[startAddressOfCorePixel + 3] = (pixels[startAddressOfPixel + 3] - depth_min) / (depth_max - depth_min) * 0.5f + 0.5f;	// range [0.5, 1]
					convert_depth_max = (convert_depth_max < corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_max;
					convert_depth_min = (convert_depth_min > corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_min;
				}
			}
		}

		// add cylinder texture
		for (int y = ymax - 9; y <= ymax; y++)
		{
			for (int x = xmin; x <= xmax; x++)
			{
				int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
				float tmp = (float(x - xmin) / (xmax - xmin) - 0.5f) * 2.0f;	//	[-1, 1]
				float t = float(x - xmin) / (xmax - xmin); // [0, 1]
				corePixels[startAddressOfCorePixel + 0] = 0.5f;
				/*corePixels[startAddressOfCorePixel + 1] = (-sqrt(1.0f - 0.25f * tmp * tmp) + 1.0f) * 0.5f;
				corePixels[startAddressOfCorePixel + 2] = (-0.5f * tmp + 1.0f) * 0.5f;*/

				corePixels[startAddressOfCorePixel + 1] = (tmp + 1.0f) * 0.5f;
				corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - tmp * tmp) + 1.0f) * 0.5f;

				//corePixels[startAddressOfCorePixel + 1] = (-0.5f + t + 1.0f) * 0.5f;
				//corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - (0.5f - t)*(0.5f - t)) + 1.0f) * 0.5f; 

				corePixels[startAddressOfCorePixel + 3] = sqrt(1.0f - tmp * tmp) * 0.5f + 0.5f;
			}
		}

		std::vector<float> textureData;
		{
			textureData.insert(textureData.end(), corePixels.begin(), corePixels.end());
		}

		glGenTextures(1, &CoreAOTexs[num_angle]);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[num_angle]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_core_AO_texture_width, g_core_AO_texture_height, 0, GL_RGBA, GL_FLOAT, textureData.data());
		glGenerateMipmap(GL_TEXTURE_2D);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		num_angle++;

		images.push_back(std::move(corePixels));
	}

	// make sure we do load some png images
	if (images.empty()) {
		std::cerr << "No images loaded." << std::endl;
		return;
	}

	std::vector<float> textureData;
	for (const auto& image : images) {
		textureData.insert(textureData.end(), image.begin(), image.end());
	}

	glGenTextures(1, &CoreAOTex);
	glBindTexture(GL_TEXTURE_3D, CoreAOTex);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_LUMINANCE, g_core_AO_texture_width, g_core_AO_texture_height, images.size(), 0, GL_LUMINANCE, GL_FLOAT, textureData.data());
	glGenerateMipmap(GL_TEXTURE_3D);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);

}

void GLWidget::InitCoreAOTexture_Long() {
	std::cout << "Read CoreAO texture: ";

	std::vector<std::string> filenames = {
	   "./high_res_AO/AO_1_Long.png",
	   "./high_res_AO/AO_1_Long.png",
	   "./high_res_AO/AO_2_Long.png",
	   "./high_res_AO/AO_3_Long.png",
	   "./high_res_AO/AO_4_Long.png",
	   "./high_res_AO/AO_5_Long.png",
	   "./high_res_AO/AO_6_Long.png",
	   "./high_res_AO/AO_7_Long.png",
	   "./high_res_AO/AO_7_Long.png"
	};
	std::vector<std::vector<float>> images;
	unsigned width = 0, height = 0;
	unsigned xmin_R, ymin_R;

	int num_angle = 0;

	std::vector<unsigned char> image_z;
	unsigned img_width_z = 0, img_height_z = 0;
	std::string z_filename = "D:/StitchMesh4.0/high_res_AO/offset_z.png";
	if (!loadRGBPNG(z_filename, image_z, img_width_z, img_height_z)) {
		std::cerr << "Failed to load image: " << z_filename << std::endl;
		return;
	}

	for (const auto& filename : filenames) {
		std::vector<unsigned char> image;
		unsigned img_width = 0, img_height = 0;

		if (!loadRGBPNG(filename, image, img_width, img_height)) {
			std::cerr << "Failed to load image: " << filename << std::endl;
			return;
		}

		GLfloat* pixels = new GLfloat[img_width * img_height * 4];
		for (size_t i = 0; i < image.size(); ++i) {
			pixels[i / 3 * 4 + i % 3] = gammaCorrection(image[i] / 255.0f, 1.f / 1);
			if (i % 3 == 2) pixels[i / 3 * 4 + 3] = gammaCorrection(image_z[i - 2] / 255.0f, 1.f / 1);
		}


		//get the empty pixel out of origin image
		// find min max
		int xmin = img_width, ymin = img_height;
		int xmax = 0, ymax = 0;
		float depth_max = 0.0;
		float depth_min = 255.0;
		for (int y = 0; y < img_height; y++)
			for (int x = 0; x < img_width; x++)
			{
				int startAddressOfPixel = ((y * img_width) + x) * 4;
				if ((pixels[startAddressOfPixel + 3] > 0.001))
				{
					if (xmin > x) xmin = x;
					if (ymin > y) ymin = y;
					if (xmax < x) xmax = x;
					if (ymax < y) ymax = y;
					depth_max = (depth_max < pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_max;
					depth_min = (depth_min > pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_min;
					//printf("pixel at 
				}
			}

		std::cout << "Core AO texture bounding box: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
		std::cout << "AO ILLUMINANCE RANGE: " << depth_max << " " << depth_min << std::endl;

		// clamp the border
		ymax += 10;

		std::cout << "height: " << ymax - ymin + 1 << " width: " << xmax - xmin + 1 << std::endl;


		// create core texture
		g_core_AO_texture_width = xmax - xmin + 1;
		g_core_AO_texture_height = ymax - ymin + 1;
		float convert_depth_max = 0.0;
		float convert_depth_min = 255.0;
		std::vector<GLfloat> corePixels(g_core_AO_texture_width * g_core_AO_texture_height * 4);
		float r_min = 255.0;
		float r_max = -255.0;
		float g_min = 255.0;
		float g_max = -255.0;
		float b_min = 255.0;
		float b_max = -255.0;
		for (int y = ymin; y <= ymax - 10; y++)
		{
			for (int x = xmin; x <= xmax; x++)
			{
				int startAddressOfPixel = ((y * img_width) + x) * 4;
				int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
				corePixels[startAddressOfCorePixel + 0] = pixels[startAddressOfPixel + 0];
				corePixels[startAddressOfCorePixel + 1] = pixels[startAddressOfPixel + 1];
				corePixels[startAddressOfCorePixel + 2] = pixels[startAddressOfPixel + 2];
				if (pixels[startAddressOfPixel + 3] < 0.001)
				{
					corePixels[startAddressOfCorePixel + 3] = 0.0;		// black part
				}
				else
				{
					r_min = (r_min > pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_min;
					g_min = (g_min > pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_min;
					b_min = (b_min > pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_min;

					r_max = (r_min < pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_max;
					g_max = (g_min < pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_max;
					b_max = (b_min < pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_max;

					corePixels[startAddressOfCorePixel + 3] = (pixels[startAddressOfPixel + 3] - depth_min) / (depth_max - depth_min) * 0.5f + 0.5f;	// range [0.5, 1]
					convert_depth_max = (convert_depth_max < corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_max;
					convert_depth_min = (convert_depth_min > corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_min;
				}
			}
		}

		// add cylinder texture
		for (int y = ymax - 9; y <= ymax; y++)
		{
			for (int x = xmin; x <= xmax; x++)
			{
				int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
				float tmp = (float(x - xmin) / (xmax - xmin) - 0.5f) * 2.0f;	//	[-1, 1]
				float t = float(x - xmin) / (xmax - xmin); // [0, 1]
				corePixels[startAddressOfCorePixel + 0] = 0.5f;
				/*corePixels[startAddressOfCorePixel + 1] = (-sqrt(1.0f - 0.25f * tmp * tmp) + 1.0f) * 0.5f;
				corePixels[startAddressOfCorePixel + 2] = (-0.5f * tmp + 1.0f) * 0.5f;*/

				corePixels[startAddressOfCorePixel + 1] = (tmp + 1.0f) * 0.5f;
				corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - tmp * tmp) + 1.0f) * 0.5f;

				//corePixels[startAddressOfCorePixel + 1] = (-0.5f + t + 1.0f) * 0.5f;
				//corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - (0.5f - t)*(0.5f - t)) + 1.0f) * 0.5f; 

				corePixels[startAddressOfCorePixel + 3] = sqrt(1.0f - tmp * tmp) * 0.5f + 0.5f;
			}
		}

		std::vector<float> textureData;
		{
			textureData.insert(textureData.end(), corePixels.begin(), corePixels.end());
		}

		glGenTextures(1, &CoreAOTexs_Long[num_angle]);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[num_angle]);
		//glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, g_core_AO_texture_width, g_core_AO_texture_height, 0, GL_LUMINANCE, GL_FLOAT, textureData.data());
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_core_AO_texture_width, g_core_AO_texture_height, 0, GL_RGBA, GL_FLOAT, textureData.data());
		glGenerateMipmap(GL_TEXTURE_2D);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		num_angle++;

		images.push_back(std::move(corePixels));
	}


}

void GLWidget::LoadCoreImage()
{
	std::vector<unsigned char> image;
	unsigned img_width = 0, img_height = 0;

	std::string filename = "./high_res_AO/normal_post.png";
	if (!loadRGBPNG(filename, image, img_width, img_height)) {
		std::cerr << "Failed to load image: " << filename << std::endl;
		return;
	}

	std::vector<unsigned char> image_z;
	unsigned img_width_z = 0, img_height_z = 0;
	std::string z_filename = "./high_res_AO/offset_z.png";
	if (!loadRGBPNG(z_filename, image_z, img_width_z, img_height_z)) {
		std::cerr << "Failed to load image: " << z_filename << std::endl;
		return;
	}

	GLfloat* pixels = new GLfloat[img_width * img_height * 4];
	for (size_t i = 0; i < image.size(); ++i) {
		pixels[i / 3 * 4 + i % 3] = gammaCorrection(image[i] / 255.0f, 1.f / 2.18);
		if (i % 3 == 2) pixels[i / 3 * 4 + 3] = gammaCorrection(image_z[i - 2] / 255.0f, 1.f / 2.18);
	}

	// find min max
	int xmin = img_width, ymin = img_height;
	int xmax = 0, ymax = 0;
	float depth_max = 0.0;
	float depth_min = 255.0;
	for (int y = 0; y < img_height; y++)
		for (int x = 0; x < img_width; x++)
		{
			int startAddressOfPixel = ((y * img_width) + x) * 4;
			if ((pixels[startAddressOfPixel + 3] > 0.001))
			{
				if (xmin > x) xmin = x;
				if (ymin > y) ymin = y;
				if (xmax < x) xmax = x;
				if (ymax < y) ymax = y;
				depth_max = (depth_max < pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_max;
				depth_min = (depth_min > pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_min;
				//printf("pixel at %d %d has color r=%f g=%f b=%f a=%f\n", x, y, pixels[startAddressOfPixel],pixels[startAddressOfPixel + 1], pixels[startAddressOfPixel + 2], pixels[startAddressOfPixel + 3]);
			}
		}

	std::cout << "Core normal texture bounding box: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
	std::cout << "depth: " << depth_max << " " << depth_min << std::endl;

	// clamp the border
	ymax += 10;

	std::cout << "height: " << ymax - ymin + 1 << " width: " << xmax - xmin + 1 << std::endl;

	// create core texture
	g_core_texture_width = xmax - xmin + 1;
	g_core_texture_height = ymax - ymin + 1;
	float convert_depth_max = 0.0;
	float convert_depth_min = 255.0;
	GLfloat* corePixels = new GLfloat[g_core_texture_width * g_core_texture_height * 4];
	float r_min = 255.0;
	float r_max = -255.0;
	float g_min = 255.0;
	float g_max = -255.0;
	float b_min = 255.0;
	float b_max = -255.0;
	for (int y = ymin; y <= ymax - 10; y++)
	{
		for (int x = xmin; x <= xmax; x++)
		{
			int startAddressOfPixel = ((y * img_width) + x) * 4;
			int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
			corePixels[startAddressOfCorePixel + 0] = pixels[startAddressOfPixel + 0];
			corePixels[startAddressOfCorePixel + 1] = pixels[startAddressOfPixel + 1];
			corePixels[startAddressOfCorePixel + 2] = pixels[startAddressOfPixel + 2];
			if (pixels[startAddressOfPixel + 3] < 0.001)
			{
				corePixels[startAddressOfCorePixel + 3] = 0.0;		// black part
			}
			else
			{
				r_min = (r_min > pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_min;
				g_min = (g_min > pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_min;
				b_min = (b_min > pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_min;

				r_max = (r_min < pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_max;
				g_max = (g_min < pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_max;
				b_max = (b_min < pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_max;

				corePixels[startAddressOfCorePixel + 3] = (pixels[startAddressOfPixel + 3] - depth_min) / (depth_max - depth_min) * 0.5f + 0.5f;	// range [0.5, 1]
				convert_depth_max = (convert_depth_max < corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_max;
				convert_depth_min = (convert_depth_min > corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_min;
			}
		}
	}

	// add cylinder texture
	for (int y = ymax - 9; y <= ymax; y++)
	{
		for (int x = xmin; x <= xmax; x++)
		{
			int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
			float tmp = (float(x - xmin) / (xmax - xmin) - 0.5f) * 2.0f;	//	[-1, 1]
			float t = float(x - xmin) / (xmax - xmin); // [0, 1]
			corePixels[startAddressOfCorePixel + 0] = 0.5f;
			/*corePixels[startAddressOfCorePixel + 1] = (-sqrt(1.0f - 0.25f * tmp * tmp) + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (-0.5f * tmp + 1.0f) * 0.5f;*/

			corePixels[startAddressOfCorePixel + 1] = (tmp + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - tmp * tmp) + 1.0f) * 0.5f;

			//corePixels[startAddressOfCorePixel + 1] = (-0.5f + t + 1.0f) * 0.5f;
			//corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - (0.5f - t)*(0.5f - t)) + 1.0f) * 0.5f; 

			corePixels[startAddressOfCorePixel + 3] = sqrt(1.0f - tmp * tmp) * 0.5f + 0.5f;
		}
	}

	glBindTexture(GL_TEXTURE_2D, coreTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_core_texture_width, g_core_texture_height, 0, GL_RGBA, GL_FLOAT, corePixels);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] corePixels;
	delete[] pixels;
}

void GLWidget::LoadCoreDirImage()
{
	std::vector<unsigned char> image;
	unsigned img_width = 0, img_height = 0;

	std::string filename = "./high_res_AO/tangent_post.png";
	if (!loadRGBPNG(filename, image, img_width, img_height)) {
		std::cerr << "Failed to load image: " << filename << std::endl;
		return;
	}

	std::vector<unsigned char> image_z;
	unsigned img_width_z = 0, img_height_z = 0;
	std::string z_filename = "./high_res_AO/offset_z.png";
	if (!loadRGBPNG(z_filename, image_z, img_width_z, img_height_z)) {
		std::cerr << "Failed to load image: " << z_filename << std::endl;
		return;
	}

	GLfloat* pixels = new GLfloat[img_width * img_height * 4];
	for (size_t i = 0; i < image.size(); ++i) {
		pixels[i / 3 * 4 + i % 3] = gammaCorrection(image[i] / 255.0f, 1.f / 2.18);
		if (i % 3 == 2) pixels[i / 3 * 4 + 3] = gammaCorrection(image_z[i - 2] / 255.0f, 1.f / 2.18);
	}

	// find min max
	int xmin = img_width, ymin = img_height;
	int xmax = 0, ymax = 0;
	float depth_max = 0.0;
	float depth_min = 255.0;
	for (int y = 0; y < img_height; y++)
		for (int x = 0; x < img_width; x++)
		{
			int startAddressOfPixel = ((y * img_width) + x) * 4;
			if ((pixels[startAddressOfPixel + 3] > 0.001))
			{
				if (xmin > x) xmin = x;
				if (ymin > y) ymin = y;
				if (xmax < x) xmax = x;
				if (ymax < y) ymax = y;
				depth_max = (depth_max < pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_max;
				depth_min = (depth_min > pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_min;
				//printf("pixel at %d %d has color r=%f g=%f b=%f a=%f\n", x, y, pixels[startAddressOfPixel],pixels[startAddressOfPixel + 1], pixels[startAddressOfPixel + 2], pixels[startAddressOfPixel + 3]);
			}
		}

	std::cout << "Core tangent texture bounding box: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
	std::cout << "depth: " << depth_max << " " << depth_min << std::endl;

	// clamp the border
	ymax += 10;

	std::cout << "height: " << ymax - ymin + 1 << " width: " << xmax - xmin + 1 << std::endl;

	// create core texture
	g_core_texture_width = xmax - xmin + 1;
	g_core_texture_height = ymax - ymin + 1;
	float convert_depth_max = 0.0;
	float convert_depth_min = 255.0;
	GLfloat* corePixels = new GLfloat[g_core_texture_width * g_core_texture_height * 4];
	float r_min = 255.0;
	float r_max = -255.0;
	float g_min = 255.0;
	float g_max = -255.0;
	float b_min = 255.0;
	float b_max = -255.0;
	for (int y = ymin; y <= ymax - 10; y++)
	{
		for (int x = xmin; x <= xmax; x++)
		{
			int startAddressOfPixel = ((y * img_width) + x) * 4;
			int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
			corePixels[startAddressOfCorePixel + 0] = pixels[startAddressOfPixel + 0];
			corePixels[startAddressOfCorePixel + 1] = pixels[startAddressOfPixel + 1];
			corePixels[startAddressOfCorePixel + 2] = pixels[startAddressOfPixel + 2];
			if (pixels[startAddressOfPixel + 3] < 0.001)
			{
				corePixels[startAddressOfCorePixel + 3] = 0.0;		// black part
			}
			else
			{
				r_min = (r_min > pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_min;
				g_min = (g_min > pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_min;
				b_min = (b_min > pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_min;

				r_max = (r_min < pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_max;
				g_max = (g_min < pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_max;
				b_max = (b_min < pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_max;

				corePixels[startAddressOfCorePixel + 3] = (pixels[startAddressOfPixel + 3] - depth_min) / (depth_max - depth_min) * 0.5f + 0.5f;	// range [0.5, 1]
				convert_depth_max = (convert_depth_max < corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_max;
				convert_depth_min = (convert_depth_min > corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_min;
			}
		}
	}

	// add cylinder texture
	for (int y = ymax - 9; y <= ymax; y++)
	{
		for (int x = xmin; x <= xmax; x++)
		{
			int startAddressOfCorePixel = (((y - ymin) * g_core_texture_width) + (x - xmin)) * 4;
			float tmp = (float(x - xmin) / (xmax - xmin) - 0.5f) * 2.0f;	//	[-1, 1]
			float t = float(x - xmin) / (xmax - xmin); // [0, 1]
			corePixels[startAddressOfCorePixel + 0] = 0.5f;
			/*corePixels[startAddressOfCorePixel + 1] = (-sqrt(1.0f - 0.25f * tmp * tmp) + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (-0.5f * tmp + 1.0f) * 0.5f;*/

			corePixels[startAddressOfCorePixel + 1] = (tmp + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - tmp * tmp) + 1.0f) * 0.5f;

			//corePixels[startAddressOfCorePixel + 1] = (-0.5f + t + 1.0f) * 0.5f;
			//corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - (0.5f - t)*(0.5f - t)) + 1.0f) * 0.5f; 

			corePixels[startAddressOfCorePixel + 3] = sqrt(1.0f - tmp * tmp) * 0.5f + 0.5f;
		}
	}

	glBindTexture(GL_TEXTURE_2D, fbrDirTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_core_texture_width, g_core_texture_height, 0, GL_RGBA, GL_FLOAT, corePixels);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] corePixels;
	delete[] pixels;
}

void GLWidget::UpdateCrossSectionTexture()
{
	std::cout << "Loading cross-section texture... ";
	glGenTextures(1, &crossSectionTex);

	int num = 60;
	std::vector<cyPoint4f>	crossSectionTextureData;
	crossSectionTextureData.resize(num);

	// get cross section points from offline plys 
	std::vector<std::pair<float, float>>& crossSectionSamples = offlinePlys->crossSectionSamples;

	for (int i = 0; i < num; i++)
	{
		std::pair<float, float>& p = crossSectionSamples[i];
		crossSectionTextureData[i].Set(p.first, p.second, ks::rand01(), ks::rand01());
	}

	glBindTexture(GL_TEXTURE_1D, crossSectionTex);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA32F, num, 0, GL_RGBA, GL_FLOAT, crossSectionTextureData.data());
	glGenerateMipmap(GL_TEXTURE_1D);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	std::cout << "Done\n";
}

void GLWidget::ReadFrameBuffer()
{
	const size_t texture_size = DEFAULT_WIDTH * DEFAULT_HEIGHT * 4 * 4 * 4;

	GLfloat* pixels = new GLfloat[texture_size];

	glBindTexture(GL_TEXTURE_2D, offlineCoreRenderBuffer->BufferTexID());

	// load texture to pixels
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, pixels);

	glBindTexture(GL_TEXTURE_2D, 0);

	// find min max
	int xmin = DEFAULT_WIDTH * 4, ymin = DEFAULT_HEIGHT * 4;
	int xmax = 0, ymax = 0;
	float depth_max = 0.0;
	float depth_min = 255.0;
	for (int y = 0; y < DEFAULT_HEIGHT * 4; y++)
		for (int x = 0; x < DEFAULT_WIDTH * 4; x++)
		{
			int startAddressOfPixel = ((y * DEFAULT_WIDTH * 4) + x) * 4;
			if ((pixels[startAddressOfPixel + 3] > 0.01))
			{
				if (xmin > x) xmin = x;
				if (ymin > y) ymin = y;
				if (xmax < x) xmax = x;
				if (ymax < y) ymax = y;
				depth_max = (depth_max < pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_max;
				depth_min = (depth_min > pixels[startAddressOfPixel + 3]) ? pixels[startAddressOfPixel + 3] : depth_min;
				//printf("pixel at %d %d has color r=%f g=%f b=%f a=%f\n", x, y, pixels[startAddressOfPixel],pixels[startAddressOfPixel + 1], pixels[startAddressOfPixel + 2], pixels[startAddressOfPixel + 3]);
			}
		}

	//int diff = ymax - DEFAULT_HEIGHT * 2;
	//diff = max(DEFAULT_HEIGHT*2 - ymin, diff);
	//xmin = DEFAULT_HEIGHT*2 - diff;
	//xmax = DEFAULT_HEIGHT * 2 + diff;

	std::cout << "Flyaway texture bounding box: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
	std::cout << "depth: " << depth_max << " " << depth_min << std::endl;

	// clamp the border
	xmin += 10;
	xmax -= 10;
	ymax += 10;

	std::cout << "height: " << ymax - ymin << " width: " << xmax - xmin << std::endl;

	// create core texture
	g_flyaway_texture_width = xmax - xmin;
	g_flyaway_texture_height = ymax - ymin;
	float convert_depth_max = 0.0;
	float convert_depth_min = 255.0;
	GLfloat* corePixels = new GLfloat[g_flyaway_texture_width * g_flyaway_texture_height * 4];
	float r_min = 255.0;
	float r_max = -255.0;
	float g_min = 255.0;
	float g_max = -255.0;
	float b_min = 255.0;
	float b_max = -255.0;
	for (int y = ymin; y < ymax - 10; y++)
	{
		for (int x = xmin; x < xmax; x++)
		{
			int startAddressOfPixel = ((y * DEFAULT_WIDTH * 4) + x) * 4;
			int startAddressOfCorePixel = (((y - ymin) * g_flyaway_texture_width) + (x - xmin)) * 4;
			corePixels[startAddressOfCorePixel + 0] = pixels[startAddressOfPixel + 0];
			corePixels[startAddressOfCorePixel + 1] = pixels[startAddressOfPixel + 1];
			corePixels[startAddressOfCorePixel + 2] = pixels[startAddressOfPixel + 2];
			if (pixels[startAddressOfPixel + 3] < 0.01)
			{
				corePixels[startAddressOfCorePixel + 3] = 0.0;		// black part
			}
			else
			{
				r_min = (r_min > pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_min;
				g_min = (g_min > pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_min;
				b_min = (b_min > pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_min;

				r_max = (r_min < pixels[startAddressOfPixel + 0]) ? pixels[startAddressOfPixel + 0] : r_max;
				g_max = (g_min < pixels[startAddressOfPixel + 1]) ? pixels[startAddressOfPixel + 1] : g_max;
				b_max = (b_min < pixels[startAddressOfPixel + 2]) ? pixels[startAddressOfPixel + 2] : b_max;

				corePixels[startAddressOfCorePixel + 3] = pixels[startAddressOfPixel + 3];	// range [0.5, 1]
				convert_depth_max = (convert_depth_max < corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_max;
				convert_depth_min = (convert_depth_min > corePixels[startAddressOfCorePixel + 3]) ? corePixels[startAddressOfCorePixel + 3] : convert_depth_min;
			}
		}
	}

	// add cylinder texture
	for (int y = ymax - 10; y < ymax; y++)
	{
		for (int x = xmin; x < xmax; x++)
		{
			int startAddressOfCorePixel = (((y - ymin) * g_flyaway_texture_width) + (x - xmin)) * 4;
			float tmp = (float(x - xmin) / (xmax - xmin) - 0.5f) * 2.0f;	//	[-1, 1]
			float t = float(x - xmin) / (xmax - xmin); // [0, 1]
			corePixels[startAddressOfCorePixel + 0] = 0.5f;
			/*corePixels[startAddressOfCorePixel + 1] = (-sqrt(1.0f - 0.25f * tmp * tmp) + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (-0.5f * tmp + 1.0f) * 0.5f;*/

			corePixels[startAddressOfCorePixel + 1] = (tmp + 1.0f) * 0.5f;
			corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - tmp * tmp) + 1.0f) * 0.5f;

			//corePixels[startAddressOfCorePixel + 1] = (-0.5f + t + 1.0f) * 0.5f;
			//corePixels[startAddressOfCorePixel + 2] = (sqrt(1.0f - (0.5f - t)*(0.5f - t)) + 1.0f) * 0.5f; 

			corePixels[startAddressOfCorePixel + 3] = sqrt(1.0f - tmp * tmp) * 0.5f + 0.5f;
			corePixels[startAddressOfCorePixel + 3] = 0.f;
		}
	}

	glBindTexture(GL_TEXTURE_2D, flyawayTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_flyaway_texture_width, g_flyaway_texture_height, 0, GL_RGBA, GL_FLOAT, corePixels);
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] corePixels;
	delete[] pixels;
}

void GLWidget::UpdateOfflineYarn()
{
	offlineFlyaway = new Fiber::Yarn();
	UpdateConfig2OfflineYarn(offlineFlyaway);
	offlineFlyaway->clock_wise = offlineFlyaway->clock_wise;
	offlineFlyaway->fiber_radius = 0.158f;
	///offlineFlyaway->z_step_num = 38;
	offlineFlyaway->z_step_num = 39;
	offlineFlyaway->aabb_micro_ct.pMin.z() = -0.5f * offlineFlyaway->z_step_num * offlineFlyaway->z_step_size;
	offlineFlyaway->yarn_radius = 0.02866;
	offlineFlyaway->aabb_micro_ct.pMax.z() = (0.5f * offlineFlyaway->z_step_num - 1) * offlineFlyaway->z_step_size;
	offlineFlyaway->simulate_fly_away();


	//---------------------------------------------      update depth based on regular fiber ---------------------------------------------   

	//glBindFramebuffer(GL_FRAMEBUFFER, offlineCoreRenderBuffer->buffer);
	//offlineCoreRenderBuffer->Bind();
	glBindFramebuffer(GL_FRAMEBUFFER, offlineCoreRenderBuffer->buffer);

	GLint vp[4]; glGetIntegerv(GL_VIEWPORT, vp);

	int w = vp[2];
	int h = vp[3];

	glViewport(0, 0, DEFAULT_WIDTH * 4, DEFAULT_HEIGHT * 4);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	GLfloat clearColor[] = { 0.0f, 0.0f, 0.0f, 0.0f };
	glClearBufferfv(GL_COLOR, 0, clearColor);

	glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);

	shaderOfflineCore.BindProgram();
	shaderOfflineCore.SetParamMatrix4(SHADER_PARAM_view_matrix, mI4.data());
	shaderOfflineCore.SetParamMatrix4(SHADER_PARAM_model_matrix, mIdentity.data());
	float fiber_thick = 0.00848868862;
	shaderOfflineCore.SetParam(SHADER_PARAM_fiber_thickness, fiber_thick);
	g_center_offset = 0.f;
	shaderOfflineCore.SetParam(SHADER_PARAM_center_offset, g_center_offset);

	{
		for (int i = 0; i < (int)offlineYarnVertices.size(); i++)
		{
			glBindBuffer(GL_ARRAY_BUFFER, offlineYarnBuffer[i]);
			glVertexAttribPointer(ATTRIB_VERTEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(ATTRIB_VERTEX);
			glDrawArrays(GL_LINE_STRIP, 0, (int)offlineYarnVertices[i].size());
			glDisableVertexAttribArray(ATTRIB_VERTEX);
		}
	}

	//---------------------------------------------      update based on flyaway fiber ---------------------------------------------   

	// clear offlinePlyVertices
	for (int i = 0; i > (int)offlineCoreVertices.size(); i++)
		offlineCoreVertices[i].clear();
	offlineCoreVertices.clear();
	for (int i = 0; i < (int)offlineCoreVertices.size(); i++)
		glDeleteBuffers(1, &offlineCoreBuffer[i]);
	offlineCoreBuffer.clear();

	for (int ply = 0; ply < offlineFlyaway->plys.size(); ply++) {
		int num_hair_fibers = offlineFlyaway->plys[ply].num_hair_fibers;
		int fiber_num = (int)offlineFlyaway->plys[ply].fibers.size();
		for (int f = fiber_num - num_hair_fibers; f < fiber_num; f++) {
			Fiber::Fiber& fiber = offlineFlyaway->plys[ply].fibers[f];
			std::vector<ks::vec3> fiberVertices;
			int fiber_vertex_num = (int)fiber.vertices.size();
			for (int v = 0; v < fiber_vertex_num; v++) {
				fiberVertices.push_back(ks::vec3(fiber.vertices[v].z(), fiber.vertices[v].y(), fiber.vertices[v].x()));
			}

			offlineCoreVertices.push_back(fiberVertices);
		}
	}

	offlineCoreBuffer.resize((int)offlineCoreVertices.size());
	// put offline ply data to buffer, used for core fiber texture generation
	for (int i = 0; i < (int)offlineCoreVertices.size(); i++)
	{
		glGenBuffers(1, &offlineCoreBuffer[i]);
		glBindBuffer(GL_ARRAY_BUFFER, offlineCoreBuffer[i]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(ks::vec3) * (offlineCoreVertices[i].size()), &offlineCoreVertices[i][0], GL_DYNAMIC_DRAW);
	}


	shaderOfflineCoreT.BindProgram();
	shaderOfflineCoreT.SetParamMatrix4(SHADER_PARAM_view_matrix, mI4.data());
	shaderOfflineCoreT.SetParamMatrix4(SHADER_PARAM_model_matrix, mIdentity.data());
	shaderOfflineCoreT.SetParam(SHADER_PARAM_fiber_thickness, fiber_thick);
	g_center_offset = 0.f;
	shaderOfflineCoreT.SetParam(SHADER_PARAM_center_offset, g_center_offset);

	for (int i = 0; i < (int)offlineCoreVertices.size(); i++)

	{
		glBindBuffer(GL_ARRAY_BUFFER, offlineCoreBuffer[i]);
		glVertexAttribPointer(ATTRIB_VERTEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(ATTRIB_VERTEX);
		glDrawArrays(GL_LINE_STRIP, 0, (int)offlineCoreVertices[i].size());
		glDisableVertexAttribArray(ATTRIB_VERTEX);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	std::string filename = "D:/test/flyaway.png";
	saveTexture2DAsImage(offlineCoreRenderBuffer->BufferTexID(), filename.data());

	ReadFrameBuffer();

	//offlineCoreRenderBuffer->Unbind();

	glViewport(0, 0, w, h);
}

bool LoadObj(const std::string& filePath, ObjData& objData) {
	std::ifstream file(filePath);
	if (!file.is_open()) {
		std::cerr << "Error: Unable to open file " << filePath << std::endl;
		return false;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string prefix;
		iss >> prefix;

		if (prefix == "v") { // Vertex position
			ks::vec3 vertex;
			iss >> vertex.x() >> vertex.y() >> vertex.z();
			objData.vertices.push_back(vertex);
		}
		else if (prefix == "vt") { // Texture coordinate
			ks::vec2 texCoord;
			iss >> texCoord.x() >> texCoord.y();
			objData.texCoords.push_back(texCoord);
		}
		else if (prefix == "vn") { // Vertex normal
			ks::vec3 normal;
			iss >> normal.x() >> normal.y() >> normal.z();
			objData.normals.push_back(normal);
		}
		else if (prefix == "f") { // Face
			Face face;
			for (int i = 0; i < 3; ++i) {
				std::string vertexData;
				iss >> vertexData;

				std::replace(vertexData.begin(), vertexData.end(), '/', ' ');
				std::istringstream vertexStream(vertexData);
				vertexStream >> face.vertexIndices[i] >> face.normalIndices[i];

				// Convert 1-based indices to 0-based indices
				face.vertexIndices[i]--;
				face.normalIndices[i]--;
			}
			objData.faces.push_back(face);
		}
	}

	file.close();
	return true;
}
void GLWidget::CreateIndexBuffer(const ObjData& objData) {
	std::vector<GLuint> indices;
	for (const auto& face : objData.faces) {
		indices.push_back(face.vertexIndices[0]);
		indices.push_back(face.vertexIndices[1]);
		indices.push_back(face.vertexIndices[2]);
	}

	glGenBuffers(1, &object_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, object_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

	std::cout << "Index Buffer Object created with " << indices.size() << " indices." << std::endl;

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void GLWidget::InitObject() {
	std::string filepath = "./foot_tr.obj";
	ObjData objData;
	LoadObj(filepath, objData);

	object_tri_num = objData.faces.size();
	{
		glGenBuffers(1, &object_vao);
		glBindBuffer(GL_ARRAY_BUFFER, object_vao);
		std::vector<GLfloat> array_vertex(0);
		for (int i = 0; i < objData.vertices.size(); i++) {
			array_vertex.push_back(GLfloat(objData.vertices[i].x()) * 30);
			array_vertex.push_back(GLfloat(objData.vertices[i].y()) * 30);
			array_vertex.push_back(GLfloat(objData.vertices[i].z()) * 30);
			array_vertex.push_back(GLfloat(objData.normals[i].x()));
			array_vertex.push_back(GLfloat(objData.normals[i].y()));
			array_vertex.push_back(GLfloat(objData.normals[i].z()));
		}
		glBufferData(GL_ARRAY_BUFFER, array_vertex.size() * sizeof(GLfloat), array_vertex.data(), GL_STATIC_DRAW);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3*sizeof(GLfloat)));
		//glEnableVertexAttribArray(0);


		//glBufferData(GL_ARRAY_BUFFER, head_meshAsset.vertexCount * 6 * sizeof(GLfloat), array_vertex.data(), GL_DYNAMIC_DRAW);

		// specify position attribute
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
		glEnableVertexAttribArray(0);

		// specify major_axis attribute
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE,
			6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	CreateIndexBuffer(objData);
}

bool GLWidget::InitShaders() {
	if (!shaderCheckerBoard.BuildProgramFiles("./glsl/checkerboard.vs.glsl", "./glsl/checkerboard.fs.glsl")
		|| !shaderYarnTess.BuildProgramFiles("./glsl/yarn_tess.vert", "./glsl/yarn_tess.frag", "./glsl/yarn_tess.geom", "./glsl/yarn_tess.tcs", "./glsl/yarn_tess.tes")
		|| !shaderOfflineYarn.BuildProgramFiles("./glsl/offlineYarn.vert", "./glsl/yarn_tess.frag", "./glsl/offlineYarn.geom")
		|| !shaderOfflineCore.BuildProgramFiles("./glsl/offlineYarn.vert", "./glsl/offlineCore.frag", "./glsl/offlineCoreT.geom")
		|| !shaderOfflineCoreT.BuildProgramFiles("./glsl/offlineYarn.vert", "./glsl/offlineCoreT.frag", "./glsl/offlineCoreT.geom")
		|| !shaderFrameBuffer.BuildProgramFiles("./glsl/framebuffer.vert", "./glsl/framebuffer.frag")
		|| !shaderDepthBuffer.BuildProgramFiles("./glsl/depthbuffer.vert", "./glsl/depthbuffer.frag")
		|| !shaderYarnTubeDepth.BuildProgramFiles("./glsl/yarn_tess.vert", "./glsl/depth.frag", "./glsl/yarn_tube.geom", "./glsl/yarn_tube.tcs", "./glsl/yarn_tube.tes")
		|| !shaderYarnTube.BuildProgramFiles("./glsl/yarn_tess.vert", "./glsl/yarn_tube.frag", "./glsl/yarn_tube.geom", "./glsl/yarn_tube.tcs", "./glsl/yarn_tube.tes")
		|| !shaderFlyaway.BuildProgramFiles("./glsl/yarn_tess.vert", "./glsl/yarn_tube_flyaway.frag", "./glsl/yarn_tube_flyaway.geom", "./glsl/yarn_tube.tcs", "./glsl/yarn_tube_flyaway.tes")
		|| !shaderYarnVerts.BuildProgramFiles("./glsl/color.vert", "./glsl/hairmeshverts.frag")
		|| !shaderDownSample.BuildProgramFiles("./glsl/downsample.vert", "./glsl/downsample.frag")
		|| !shaderobject.BuildProgramFiles("./glsl/head.vert", "./glsl/head.frag")
		|| !shaderobjectshadow.BuildProgramFiles("./glsl/head.vert")
		) return false;

	shaderFlyaway.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderFlyaway.RegisterParam(SHADER_PARAM_camera_matrix, "camera_matrix");
	shaderFlyaway.RegisterParam(SHADER_PARAM_shadow_matrix, "shadow_matrix");
	shaderFlyaway.RegisterParam(SHADER_PARAM_shadow_texture, "shadow_tex");
	shaderFlyaway.RegisterParam(SHADER_PARAM_flyaway_texture, "flyaway_tex");
	shaderFlyaway.RegisterParam(SHADER_PARAM_light_dir, "light_dir");
	shaderFlyaway.RegisterParam(SHADER_PARAM_view_dir, "view_dir");
	shaderFlyaway.RegisterParam(SHADER_PARAM_AR, "AR");
	shaderFlyaway.RegisterParam(SHADER_PARAM_ATT, "ATT");
	shaderFlyaway.RegisterParam(SHADER_PARAM_AD, "AD");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Beta_R, "betaR");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Beta_TT, "betaTT");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Beta_N, "betaN");
	shaderFlyaway.RegisterParam(SHADER_PARAM_D_B, "d_b");
	shaderFlyaway.RegisterParam(SHADER_PARAM_D_F, "d_f");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Beta_Diffuse, "beta_diffuse");
	shaderFlyaway.RegisterParam(SHADER_PARAM_D_F_Inner, "d_f_inner");
	shaderFlyaway.RegisterParam(SHADER_PARAM_D_Cross_Inter, "d_cross_inter");
	shaderFlyaway.RegisterParam(SHADER_PARAM_TUBE_WIDTH, "tube_width");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Delta_2, "width_scale");
	shaderFlyaway.RegisterParam(SHADER_PARAM_Delta_0, "mDelta_0");
	shaderFlyaway.RegisterParam(SHADER_PARAM_flyaway_texture_height, "flyaway_texture_height");
	shaderFlyaway.RegisterParam(SHADER_PARAM_TEXTURE_OFFSET, "texture_offset");
	shaderFlyaway.BindProgram();
	shaderFlyaway.SetParam(SHADER_PARAM_shadow_texture, TEXTURE_UNIT_SHADOW_TEXTURE);
	shaderFlyaway.SetParam(SHADER_PARAM_flyaway_texture, TEXTURE_UNIT_FLYAWAY_TEXTURE);

	shaderDownSample.RegisterParam(SHADER_PARAM_SSAA_texture, "uTexture");
	shaderDownSample.RegisterParam(SHADER_PARAM_SSAA_texture2, "uTexture2");
	shaderDownSample.BindProgram();
	shaderDownSample.SetParam(SHADER_PARAM_SSAA_texture, TEXTURE_UNIT_SSAA_TEXTURE);
	shaderDownSample.SetParam(SHADER_PARAM_SSAA_texture2, TEXTURE_UNIT_SSAA_TEXTURE2);

	shaderYarnVerts.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderYarnVerts.RegisterParam(SHADER_PARAM_model_matrix, "modelMatrix");
	//shaderYarnVerts.RegisterParam(SHADER_PARAM_camera_matrix, "camera_matrix");
	//shaderYarnVerts.RegisterParam(SHADER_PARAM_light_pos_world, "light_pos_world");

	// yarn tessellation shader
	shaderYarnTess.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderYarnTess.RegisterParam(SHADER_PARAM_shadow_matrix, "shadow_matrix");
	shaderYarnTess.RegisterParam(SHADER_PARAM_camera_matrix, "camera_matrix");
	shaderYarnTess.RegisterParam(SHADER_PARAM_shadow_R_matrix, "shadow_R_matrix");

	shaderYarnTess.RegisterParam(SHADER_PARAM_light_dir, "light_dir");
	//shaderYarnTess.RegisterParam(SHADER_PARAM_light_pos, "light_pos");
	shaderYarnTess.RegisterParam(SHADER_PARAM_light_pos_world, "light_pos_world");
	shaderYarnTess.RegisterParam(SHADER_PARAM_view_dir, "view_dir");
	shaderYarnTess.RegisterParam(SHADER_PARAM_view_pos, "view_pos");

	shaderYarnTess.RegisterParam(SHADER_PARAM_color, "default_color");
	shaderYarnTess.RegisterParam(SHADER_PARAM_lamda_R, "lamda_R");
	shaderYarnTess.RegisterParam(SHADER_PARAM_k_R, "k_R");

	shaderYarnTess.RegisterParam(SHADER_PARAM_cross_section_texture, "cross_section_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_lod, "use_lod");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_lod_vis, "use_lod_vis");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_ao, "use_ao");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_diffuse, "use_diffuse");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_specular, "use_specular");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_regular_fiber, "use_regular_fiber");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_self_shadow, "use_self_shadow");
	shaderYarnTess.RegisterParam(SHADER_PARAM_use_shadow, "use_shadow");
	shaderYarnTess.RegisterParam(SHADER_PARAM_inv_ply_num, "inv_ply_num");

	// fiber parameters
	shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_thickness, "fiber_thickness");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_rho_min, "rho_min");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_rho_max, "rho_max");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_s_i, "s_i");
	shaderYarnTess.RegisterParam(SHADER_PARAM_inv_fiber_alpha, "inv_alpha_fiber");
	//shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_use_migration, "use_migration");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fiber_use_fly_away, "use_fly_away");

	// ply parameters
	shaderYarnTess.RegisterParam(SHADER_PARAM_ply_ellipse_long, "ellipse_long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_ply_ellipse_short, "ellipse_short");
	shaderYarnTess.RegisterParam(SHADER_PARAM_inv_ply_alpha, "inv_alpha_ply");
	shaderYarnTess.RegisterParam(SHADER_PARAM_ply_radius, "radius_ply");
	shaderYarnTess.RegisterParam(SHADER_PARAM_ply_use_core_fibers, "use_core_flys");

	// fly away
	shaderYarnTess.RegisterParam(SHADER_PARAM_fly_away_r_loop_max, "radius_loop_max");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fly_away_rho_loop, "rho_loop");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fly_away_rho_hair, "rho_hair");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fly_away_hair_rot_scale, "hair_rot_scale");
	shaderYarnTess.RegisterParam(SHADER_PARAM_fly_away_hair_len_scale, "hair_len_scale");
	shaderYarnTess.RegisterParam(SHADER_PARAM_hair_scale, "hair_scale");

	shaderYarnTess.RegisterParam(SHADER_PARAM_ply_use_core_texture, "use_core_texture");

	//shaderYarnTess.RegisterParam(SHADER_PARAM_scale, "scale");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_texture_height, "core_texture_height");

	shaderYarnTess.RegisterParam(SHADER_PARAM_self_shadow_texture, "self_shadow_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_shadow_texture, "shadow_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_texture, "core_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_dir_texture, "fbr_dir_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO_texture, "Core_AO_tex");

	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO0_texture, "Core_AO0_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO1_texture, "Core_AO1_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO2_texture, "Core_AO2_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO3_texture, "Core_AO3_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO4_texture, "Core_AO4_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO5_texture, "Core_AO5_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO6_texture, "Core_AO6_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO7_texture, "Core_AO7_tex");

	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO0_texture_Long, "Core_AO0_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO1_texture_Long, "Core_AO1_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO2_texture_Long, "Core_AO2_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO3_texture_Long, "Core_AO3_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO4_texture_Long, "Core_AO4_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO5_texture_Long, "Core_AO5_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO6_texture_Long, "Core_AO6_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO7_texture_Long, "Core_AO7_tex_Long");
	shaderYarnTess.RegisterParam(SHADER_PARAM_core_AO8_texture_Long, "Core_AO8_tex_Long");

	shaderYarnTess.RegisterParam(SHADER_PARAM_shadow_R_texture, "shadow_R_tex");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Delta_0, "mDelta_0");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Delta_1, "mDelta_1");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Delta_2, "mDelta_2");
	shaderYarnTess.RegisterParam(SHADER_PARAM_AR, "AR");
	shaderYarnTess.RegisterParam(SHADER_PARAM_ATT, "ATT");
	shaderYarnTess.RegisterParam(SHADER_PARAM_AD, "AD");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Beta_R, "betaR");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Beta_TT, "betaTT");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Beta_N, "betaN");
	shaderYarnTess.RegisterParam(SHADER_PARAM_D_B, "d_b");
	shaderYarnTess.RegisterParam(SHADER_PARAM_D_F, "d_f");
	shaderYarnTess.RegisterParam(SHADER_PARAM_Beta_Diffuse, "beta_diffuse");
	shaderYarnTess.RegisterParam(SHADER_PARAM_D_F_Inner, "d_f_inner");
	shaderYarnTess.RegisterParam(SHADER_PARAM_D_Cross_Inter, "d_cross_inter");
	shaderYarnTess.RegisterParam(SHADER_PARAM_GRID_RES, "grid_res");
	shaderYarnTess.RegisterParam(SHADER_PARAM_GRID_RADIUS, "grid_radius");
	shaderYarnTess.RegisterParam(SHADER_PARAM_SCENE_MIN, "scene_min");
	shaderYarnTess.RegisterParam(SHADER_PARAM_DENSITY_TEXTURE, "density_tex");

	shaderYarnTess.BindProgram();
	shaderYarnTess.SetParam(SHADER_PARAM_cross_section_texture, TEXTURE_UNIT_CROSS_SECTION_TEXUTRE);
	shaderYarnTess.SetParam(SHADER_PARAM_self_shadow_texture, TEXTURE_UNIT_SELFSHADOW_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_shadow_texture, TEXTURE_UNIT_SHADOW_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_DENSITY_TEXTURE, TEXTURE_UNIT_DENSITY_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_texture, TEXTURE_UNIT_CORE_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_dir_texture, TEXTURE_UNIT_CORE_DIR_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO_texture, TEXTURE_UNIT_CORE_AO_TEXTURE);

	shaderYarnTess.SetParam(SHADER_PARAM_core_AO0_texture_Long, TEXTURE_UNIT_CORE_AO0_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO1_texture_Long, TEXTURE_UNIT_CORE_AO1_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO2_texture_Long, TEXTURE_UNIT_CORE_AO2_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO3_texture_Long, TEXTURE_UNIT_CORE_AO3_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO4_texture_Long, TEXTURE_UNIT_CORE_AO4_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO5_texture_Long, TEXTURE_UNIT_CORE_AO5_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO6_texture_Long, TEXTURE_UNIT_CORE_AO6_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO7_texture_Long, TEXTURE_UNIT_CORE_AO7_TEXTURE_Long);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO8_texture_Long, TEXTURE_UNIT_CORE_AO8_TEXTURE_Long);

	shaderYarnTess.SetParam(SHADER_PARAM_core_AO0_texture, TEXTURE_UNIT_CORE_AO0_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO1_texture, TEXTURE_UNIT_CORE_AO1_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO2_texture, TEXTURE_UNIT_CORE_AO2_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO3_texture, TEXTURE_UNIT_CORE_AO3_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO4_texture, TEXTURE_UNIT_CORE_AO4_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO5_texture, TEXTURE_UNIT_CORE_AO5_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO6_texture, TEXTURE_UNIT_CORE_AO6_TEXTURE);
	shaderYarnTess.SetParam(SHADER_PARAM_core_AO7_texture, TEXTURE_UNIT_CORE_AO7_TEXTURE);


	shaderYarnTess.SetParam(SHADER_PARAM_shadow_R_texture, TEXTURE_UNIT_SHADOW_R_TEXTURE);


	// checker board shader
	shaderCheckerBoard.RegisterParam(SHADER_PARAM_checker_mv_matrix, "mv_matrix");
	shaderCheckerBoard.RegisterParam(SHADER_PARAM_checker_proj_matrix, "proj_matrix");
	shaderCheckerBoard.RegisterParam(SHADER_PARAM_checker_texture, "tex_object");
	shaderCheckerBoard.RegisterParam(SHADER_PARAM_light_pos, "light_pos");
	shaderCheckerBoard.BindProgram();
	shaderCheckerBoard.SetParam(SHADER_PARAM_checker_texture, TEXTURE_UNIT_CHECKER_TEXUTRE);


	// offline yarn shader
	//shaderOfflineYarn.RegisterParam(SHADER_PARAM_light_dir, "light_dir");
	//shaderOfflineYarn.RegisterParam(SHADER_PARAM_light_pos, "lightPos");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_model_matrix, "modelMatrix");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_camera_matrix, "camera_matrix");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_light_pos_world, "light_pos_world");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_view_dir, "view_dir");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_fiber_thickness, "fiber_thickness");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_use_ao, "use_ao");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_use_diffuse, "use_diffuse");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_use_self_shadow, "use_self_shadow");
	shaderOfflineYarn.RegisterParam(SHADER_PARAM_center_offset, "center_offset");
	//shaderOfflineYarn.RegisterParam(SHADER_PARAM_scale, "scale");

	shaderOfflineYarn.RegisterParam(SHADER_PARAM_color, "default_color");

	shaderOfflineYarn.RegisterParam(SHADER_PARAM_self_shadow_texture, "self_shadow_tex");

	shaderOfflineYarn.BindProgram();
	shaderOfflineYarn.SetParam(SHADER_PARAM_self_shadow_texture, TEXTURE_UNIT_SELFSHADOW_TEXTURE);


	// offline core shader
	shaderOfflineCore.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderOfflineCore.RegisterParam(SHADER_PARAM_model_matrix, "modelMatrix");
	shaderOfflineCore.RegisterParam(SHADER_PARAM_fiber_thickness, "fiber_thickness");
	shaderOfflineCore.RegisterParam(SHADER_PARAM_center_offset, "center_offset");

	// offline core tangent shader
	shaderOfflineCoreT.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderOfflineCoreT.RegisterParam(SHADER_PARAM_model_matrix, "modelMatrix");
	shaderOfflineCoreT.RegisterParam(SHADER_PARAM_fiber_thickness, "fiber_thickness");
	shaderOfflineCoreT.RegisterParam(SHADER_PARAM_center_offset, "center_offset");


	shaderYarnTubeDepth.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderYarnTubeDepth.RegisterParam(SHADER_PARAM_TUBE_WIDTH, "tube_width");
	shaderYarnTubeDepth.RegisterParam(SHADER_PARAM_Light_Dir, "light_dir");
	shaderYarnTubeDepth.BindProgram();


	// frame buffer shader
	shaderFrameBuffer.RegisterParam(SHADER_PARAM_frame_buffer_texture, "framebuffer_tex");
	shaderFrameBuffer.BindProgram();
	shaderFrameBuffer.SetParam(SHADER_PARAM_frame_buffer_texture, TEXTURE_UNIT_FRAMEBUFFER_TEXUTRE);


	shaderDepthBuffer.RegisterParam(SHADER_PARAM_depth_buffer_texture, "depthbuffer_tex");
	shaderDepthBuffer.BindProgram();
	shaderDepthBuffer.SetParam(SHADER_PARAM_depth_buffer_texture, TEXTURE_UNIT_DEPTHBUFFER_TEXUTRE);

	shaderYarnTube.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderYarnTube.RegisterParam(SHADER_PARAM_shadow_matrix, "shadow_matrix");
	shaderYarnTube.RegisterParam(SHADER_PARAM_shadow_texture, "shadow_tex");
	shaderYarnTube.BindProgram();
	shaderYarnTube.SetParam(SHADER_PARAM_shadow_texture, TEXTURE_UNIT_SHADOW_TEXTURE);

	shaderobjectshadow.RegisterParam(SHADER_PARAM_shadow_matrix, "view_matrix");

	shaderobject.RegisterParam(SHADER_PARAM_view_matrix, "view_matrix");
	shaderobject.RegisterParam(SHADER_PARAM_shadow_matrix, "shadow_matrix");
	shaderobject.RegisterParam(SHADER_PARAM_shadow_texture, "shadow_tex");
	shaderobject.RegisterParam(SHADER_PARAM_light_dir, "light_dir");
	shaderobject.RegisterParam(SHADER_PARAM_view_pos, "cameraPos");
	shaderobject.BindProgram();
	shaderobject.SetParam(SHADER_PARAM_shadow_texture, TEXTURE_UNIT_SHADOW_TEXTURE);

	return true;
}

void glfw_error_callback(int error, const char* description) {
	fprintf(stderr, "GLFW Error (%d): %s\n", error, description);
}

void GLWidget::initializeGL()
{
	// Setup window
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) {
		printf("glfw unable to initialize");
		return;
	}

	const char* glsl_version = "#version 460";
	int SSAA = 1;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glEnable(GL_MULTISAMPLE);
	//glDisable(GL_MULTISAMPLE); 
	glfwWindowHint(GLFW_SAMPLES, SSAA);
	//glfwWindowHint(GLFW_SAMPLES, SSAA);

	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "GLFW+GLEW+IMGUI Cloth Rendering", NULL, NULL);
	if (window == NULL) {
		printf("can not create a window");
		return;
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(0); // disable vsync

	// ImGUI process
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	if (glewInit() != GLEW_OK)
	{
		std::cout << "Error" << std::endl;
	}

	glClearColor(0.2f, 0.2f, 0.2f, 0);

	glPointSize(3.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_CULL_FACE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_MULTISAMPLE);

	glGetIntegerv(GL_CURRENT_PROGRAM, &DefaultShaderProgram);

	if (!InitShaders()) assert(0);

	mIdentity = Eigen::Matrix4f::Identity();
	mI4 = mIdentity;
	mI4(0, 0) = 4;
	mI4(1, 1) = 4;
	mI4(2, 2) = 4;

	//Initialization of buffers
	InitOfflineYarns();
	UpdateCrossSectionTexture();
	Init1DCylinderNormTexture();
	InitObject();

	fiberData.g_yarn_radius = 0.02866;
	fiberData.g_fiber_thickness = 0.008f;

	offlineCoreRenderBuffer = new OpenGLRenderBuffer;
	offlineCoreRenderBuffer->Initialize(true);
	offlineCoreRenderBuffer->Resize(DEFAULT_WIDTH * 4, DEFAULT_HEIGHT * 4, OpenGLRenderBuffer::FORMAT_FLOAT_32, 4);

	offlineCoreTangentRenderBuffer = new OpenGLRenderBuffer;
	offlineCoreTangentRenderBuffer->Initialize(true);
	offlineCoreTangentRenderBuffer->Resize(DEFAULT_WIDTH * 4, DEFAULT_HEIGHT * 4, OpenGLRenderBuffer::FORMAT_FLOAT_32, 4);

	int w = WINDOW_WIDTH;
	int h = WINDOW_HEIGHT;

	SSAARenderBuffer = new OpenGLRenderBuffer;
	SSAARenderBuffer->Initialize(true, 2);
	SSAARenderBuffer->Resize(w * 2, h * 2, OpenGLRenderBuffer::FORMAT_FLOAT_32, 4);

	SSAARenderBuffer2 = new OpenGLRenderBuffer;
	SSAARenderBuffer2->Initialize(true, 2);
	SSAARenderBuffer2->Resize(w * 2, h * 2, OpenGLRenderBuffer::FORMAT_FLOAT_32, 4);

	InitShadowTextureBuffer();

	glGenTextures(1, &coreTex);
	//UpdateCoreTexture();
	LoadCoreImage();

	glGenTextures(1, &fbrDirTex);
	//UpdateCoreDirTexture();
	LoadCoreDirImage();

	glGenTextures(1, &flyawayTex);
	UpdateOfflineYarn();

	glGenTextures(1, &CoreAOTex);
	InitCoreAOTexture();

	InitCoreAOTexture_Long();
	//fiberData.g_color[0] = g_color_chart[g_shading_idx - 1][0];
	//fiberData.g_color[1] = g_color_chart[g_shading_idx - 1][1];
	//fiberData.g_color[2] = g_color_chart[g_shading_idx - 1][2];

	g_shading_idx = 1;
}

float bsplineBase3(int idx, double u)
{
	double rtn = 0;
	switch (idx) {
	case 0:
		rtn = (1 - u) * (1 - u) * (1 - u) / 6;
		break;
	case 1:
		rtn = (3 * u * u * u - 6 * u * u + 4) / 6;
		break;
	case 2:
		rtn = (1 + 3 * u + 3 * u * u - 3 * u * u * u) / 6;
		break;
	case 3:
		rtn = u * u * u / 6;
		break;
	}
	return rtn;
}

inline ks::vec3 cubicQxy(float t, ks::vec3 c0, ks::vec3 c1, ks::vec3 c2, ks::vec3 c3)
{
	return c0 * bsplineBase3(0, t) + c1 * bsplineBase3(1, t) + c2 * bsplineBase3(2, t) + c3 * bsplineBase3(3, t);
}

float cubicBezierLength(ks::vec3& c0, ks::vec3& c1, ks::vec3& c2, ks::vec3& c3)
{
	int ptCount = 40;
	float totDist = 0;
	ks::vec3 lastX = c0;
	for (int i = 1; i <= ptCount; i++) {
		ks::vec3 pt = cubicQxy(i / ptCount, c0, c1, c2, c3);
		ks::vec3 d = pt.head(3) - lastX.head(3);
		totDist += d.norm();
		lastX = pt;
	}
	return totDist;
}

void GLWidget::InitShadowTextureBuffer()
{
	int prevBuffer;
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevBuffer);

	//depth_fbo = new QOpenGLFramebufferObject(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE, QOpenGLFramebufferObject::Depth);
	glGenFramebuffers(1, &depth_fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, depth_fbo);

	glGenTextures(1, &depth_tex);
	glBindTexture(GL_TEXTURE_2D, depth_tex);

	glTexStorage2D(GL_TEXTURE_2D, 11, GL_DEPTH_COMPONENT32F, SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_tex, 0);

	glGenTextures(1, &depth_debug_tex);
	glBindTexture(GL_TEXTURE_2D, depth_debug_tex);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32F, SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, depth_debug_tex, 0);


	glBindFramebuffer(GL_FRAMEBUFFER, prevBuffer);

	glEnable(GL_DEPTH_TEST);
}

void GLWidget::UpdateShadowTexture()
{
	//specific for flame pattern
	//g_light_pos = ks::vec3(0, 0, 40);
	znear = 15.f;
	zfar = 60.f;
	fov = 60.f;


	//light_proj_matrix = ks::SetPerspective(DEG2RAD(fov), aspect, znear, zfar);

	//shadow pass orthogonal projection parameter
	float t, b, l, r;
	////for glove
	t = 8;
	b = -8;
	l = -8;
	r = 8;

	//for woven pattern
	//t = 4;
	//b = -4;
	//l = -4;
	//r = 4;
	znear = 0.5;
	zfar = 400;

	light_proj_matrix = ks::SetOrtho(l, r, t, b, znear, zfar);
	light_view_matrix = ks::SetView(g_light_pos + center, center, ks::vec3(0.0f, 1.0f, 0.0f));
	//light_view_matrix = mIdentity;

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glEnable(GL_DEPTH_TEST);

	glDisable(GL_CULL_FACE);

	shaderYarnTubeDepth.BindProgram();
	ks::mat4 light_mat = (light_proj_matrix * light_view_matrix);
	shaderYarnTubeDepth.SetParamMatrix4(SHADER_PARAM_view_matrix, light_mat.data());
	shaderYarnTubeDepth.SetParam(SHADER_PARAM_TUBE_WIDTH, fiberData.g_yarn_radius);
	ks::vec3 light_Dir = g_light_pos.normalized();
	shaderYarnTubeDepth.SetParam(SHADER_PARAM_Light_Dir, light_Dir.x(), light_Dir.y(), light_Dir.z());

	glBindBuffer(GL_ARRAY_BUFFER, vboYarn[VBO_YARN_VERTEX]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * array_Vertex.size(), NULL, GL_DYNAMIC_DRAW);

	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * array_Vertex.size(), array_Vertex.data());

	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);

	ks::mat4 mat = light_proj_matrix * light_view_matrix;
	//transpose before set into opengl
	ks::mat4 gl_mat = mat;
	shaderYarnTubeDepth.SetParamMatrix4(SHADER_PARAM_view_matrix, gl_mat.data());

	glPatchParameteri(GL_PATCH_VERTICES, 4);
	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(4.0f, 4.0f);

	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);
	//glEnable(GL_DEPTH_TEST);

	//glDisable(GL_CULL_FACE);

#if USE_KNIT
	glDrawArrays(GL_PATCHES, 0, (int)_knit[currModel][g_current_frame_num]._long_yarn_ctl_pts.size());
#else
	glDrawArrays(GL_PATCHES, 0, vboYarnVertCount);
	//glDrawArrays(GL_TRIANGLES, 0, vboYarnVertCount/3 * 3);
#endif

	glDisable(GL_POLYGON_OFFSET_FILL);

}


void GLWidget::render_task(float* Yarn_ctrPoints, int* first_ctrP_idx, int yarn_num)
{
	array_Vertex.resize(0);
	ks::AABB3 bound;
	bound.reset();
	for (int yarn_idx = 0; yarn_idx < yarn_num - 1; yarn_idx++) {
		for (int i = first_ctrP_idx[yarn_idx]; i < first_ctrP_idx[yarn_idx + 1] - 3; i += 1) {
			ks::vec3 c0 = ks::vec3(Yarn_ctrPoints[i * 3], Yarn_ctrPoints[i * 3 + 1], Yarn_ctrPoints[i * 3 + 2]);
			ks::vec3 c1 = ks::vec3(Yarn_ctrPoints[(i + 1) * 3], Yarn_ctrPoints[(i + 1) * 3 + 1], Yarn_ctrPoints[(i + 1) * 3 + 2]);
			ks::vec3 c2 = ks::vec3(Yarn_ctrPoints[(i + 2) * 3], Yarn_ctrPoints[(i + 2) * 3 + 1], Yarn_ctrPoints[(i + 2) * 3 + 2]);
			ks::vec3 c3 = ks::vec3(Yarn_ctrPoints[(i + 3) * 3], Yarn_ctrPoints[(i + 3) * 3 + 1], Yarn_ctrPoints[(i + 3) * 3 + 2]);

			ks::vec3 prevOne = ks::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
			float segmentArcLen = 0;
			for (int j = 0; j <= 10; j++) {
				const float t = float(j) / float(10);
				// cubic b-spline
				ks::vec3 p = cubicQxy(t, c0, c1, c2, c3);
				if (prevOne.x() != FLT_MAX)
					segmentArcLen += (p - prevOne).norm();
				prevOne = p;
			}
		}
		// estimation of number of segments per patch
		int segment = 10; // segment per patch

		float worldSpaceStepSize = 0.01;
		std::vector<ks::vec3> yarn_divided(0);

		float arclen = 0;
		// step through arclen
		// std::cout << "segment is " << segment << " " << newsize << " " << numCurvePatches <<
		// std::endl;
		float carryThrough = 0.f;
		for (int i = first_ctrP_idx[yarn_idx]; i < first_ctrP_idx[yarn_idx + 1] - 3; i += 1) {
			ks::vec3 c0 = ks::vec3(Yarn_ctrPoints[i * 3], Yarn_ctrPoints[i * 3 + 1], Yarn_ctrPoints[i * 3 + 2]);
			ks::vec3 c1 = ks::vec3(Yarn_ctrPoints[(i + 1) * 3], Yarn_ctrPoints[(i + 1) * 3 + 1], Yarn_ctrPoints[(i + 1) * 3 + 2]);
			ks::vec3 c2 = ks::vec3(Yarn_ctrPoints[(i + 2) * 3], Yarn_ctrPoints[(i + 2) * 3 + 1], Yarn_ctrPoints[(i + 2) * 3 + 2]);
			ks::vec3 c3 = ks::vec3(Yarn_ctrPoints[(i + 3) * 3], Yarn_ctrPoints[(i + 3) * 3 + 1], Yarn_ctrPoints[(i + 3) * 3 + 2]);

			ks::vec3 prevOne = ks::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
			float segmentArcLen = 0;
			for (int j = 0; j <= segment; j++) {
				const float t = float(j) / float(segment);
				// cubic b-spline
				ks::vec3 p = cubicQxy(t, c0, c1, c2, c3);
				if (prevOne.x() != FLT_MAX)
					segmentArcLen += (p - prevOne).norm();
				prevOne = p;
			}
			arclen += segmentArcLen;

			float nextPos = carryThrough;
			while (nextPos < segmentArcLen) {
				float t = nextPos / segmentArcLen;

				// std::cout << "nextPos is " << nextPos << "segmentarc " << segmentArcLen << std::endl;
				//  put one
				//  cubic b-spline
				ks::vec3 cur_t = cubicQxy(t, c0, c1, c2, c3);
				ks::vec3 cur_t_1 = cubicQxy(t + 0.01f, c0, c1, c2, c3);

				yarn_divided.push_back(cur_t);

				nextPos += worldSpaceStepSize;
			}
			carryThrough = nextPos - segmentArcLen;

			//if (i + 3 >= Yarn_file.firstControlPoint_[yarn_idx + 1] - 3) {
			//	// a rounding method that "on average" gets the right length
			//	if (carryThrough < 0.5f * worldSpaceStepSize) {
			//		float t = nextPos / segmentArcLen;
			//		vec3 cur_t = cubicQxy(t, c0, c1, c2, c3);
			//		vec3 cur_t_1 = cubicQxy(t + 0.01f, c0, c1, c2, c3);
			//		yarn_divided.push_back(cur_t);
			//	}
			//}
		}
		//std::cout << yarn_idx << " 's length is " << arclen / 0.3802 << std::endl;
		arclen = 0;
		for (int i = 0; i < yarn_divided.size() - 3; i++) {

			ks::vec3 c0, c1, c2, c3;
			{
				c0 = yarn_divided[i];
				c1 = yarn_divided[i + 1];
				c2 = yarn_divided[i + 2];
				c3 = yarn_divided[i + 3];
			}
			array_Vertex.push_back(c0.x());
			array_Vertex.push_back(c0.y());
			array_Vertex.push_back(c0.z());
			array_Vertex.push_back(arclen);
			bound.expand(c0);

			ks::vec3 prevOne = ks::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
			float segmentArcLen = 0;
			for (int j = 0; j <= 10; j++) {
				const float t = float(j) / float(10);
				// cubic b-spline
				ks::vec3 p = cubicQxy(t, c0, c1, c2, c3);
				if (prevOne.x() != FLT_MAX)
					segmentArcLen += (p - prevOne).norm();
				prevOne = p;
			}
			//arclen += segmentArcLen;
			//std::cout << i << "th segment length is " << segmentArcLen << std::endl;
			segmentArcLen = worldSpaceStepSize;
			segmentArcLen = segmentArcLen / 0.38;
			arclen += segmentArcLen;

			array_Vertex.push_back(c1.x());
			array_Vertex.push_back(c1.y());
			array_Vertex.push_back(c1.z());
			array_Vertex.push_back(arclen);

			array_Vertex.push_back(c2.x());
			array_Vertex.push_back(c2.y());
			array_Vertex.push_back(c2.z());
			array_Vertex.push_back(arclen);

			array_Vertex.push_back(c3.x());
			array_Vertex.push_back(c3.y());
			array_Vertex.push_back(c3.z());
			array_Vertex.push_back(arclen);

		}
		std::cout << yarn_idx << " s real arclen is " << arclen << std::endl;
	}
	vboYarnVertCount = array_Vertex.size() / 4;
	center = bound.center();
	camera.SetTarget(center);

	std::cout << "center is " << center << std::endl;
	glGenBuffers(1, &vboYarn[VBO_YARN_VERTEX]);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	ks::vec2 initialPos;
	ImVec2 last_mouse_pos = ImGui::GetMousePos();

	vgm::Quat qRot = vgm::Quat(1.f, 0.f, 0.f, 0.f);
	vgm::Vec3 lightControlDir = { 0.64f, -0.72f, 0.13f };
	float thickness = fiberData.g_yarn_radius * 10;
	while (!glfwWindowShouldClose(window))
	{

		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);

		//----------------------------------------------------   Shadow Pass ------------------------------------------------------------------------------------
		//for (int i = 0; i < lightData.lightsCount; i++) 
		{
			glViewport(0, 0, SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
			glBindFramebuffer(GL_FRAMEBUFFER, depth_fbo);
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(4.0f, 4.0f);

			static const GLenum buffs[] = { GL_COLOR_ATTACHMENT0 };
			glDrawBuffers(1, buffs);
			glClearBufferfv(GL_COLOR, 0, zero);
			glClearBufferfv(GL_DEPTH, 0, ones);

			UpdateShadowTexture();

			//shadow map about object part
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, object_ibo);
			glBindBuffer(GL_ARRAY_BUFFER, object_vao);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
			glEnableVertexAttribArray(0);

			// specify normal attribute
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE,
				6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
			glEnableVertexAttribArray(1);

			shaderobjectshadow.BindProgram();

			ks::mat4 light_mat = (light_proj_matrix * light_view_matrix);
			shaderobjectshadow.SetParamMatrix4(SHADER_PARAM_shadow_matrix, light_mat.data());

			glDrawElements(GL_TRIANGLES, object_tri_num * 3, GL_UNSIGNED_INT, NULL);
			glDisableVertexAttribArray(1);

			//yarn shadow

			glViewport(0, 0, display_w, display_h);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		//-------------------------------------------------------------- second stage----------------------------------------------------------

		int w = display_w;
		int h = display_h;
		//std::cout << " w ,h is " << w << " " << h << std::endl;
		SSAARenderBuffer->Resize(w * 3, h * 3, OpenGLRenderBuffer::FORMAT_FLOAT_32, 4);
		SSAARenderBuffer->Bind();

		glViewport(0, 0, w*3, h*3);
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
		GLfloat clearColor[] = { 0.2f, 0.2f, 0.2f, 1.f };
		glClearBufferfv(GL_COLOR, 0, clearColor);

		ks::Transform vp_transform(camera.GetMatrix());
		ks::Transform camera_trans(camera.GetMatrix());
		ks::vec3 vd = camera_trans.point_hdiv(ks::vec3(0, 0, 1));
		vd.normalized();
		ks::mat4 lookup;
		lookup << 0.5f, 0, 0, 0.5f,
			0, 0.5f, 0, 0.5f,
			0, 0, 0.5f, 0.5f,
			0, 0, 0, 1;
		ks::mat4 shadowMatrix = lookup * light_proj_matrix * light_view_matrix;
		ks::mat4 shadowMat = light_proj_matrix * light_view_matrix;

		glPatchParameteri(GL_PATCH_VERTICES, 4);

		//update camera and lighting
		ks::mat4 trans = camera.GetMatrix();
		lightPosCam4 = trans * ks::vec4(g_light_pos.x(), g_light_pos.y(), g_light_pos.z(), 1.f);

		camera_matrix_m4 = camera.GetViewMatrix();

		//texture binding
		shaderYarnTess.BindProgram();

		//glDisable(GL_MULTISAMPLE); 

		glBindBuffer(GL_ARRAY_BUFFER, vboYarn[VBO_YARN_VERTEX]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * array_Vertex.size(), NULL, GL_DYNAMIC_DRAW);

		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * array_Vertex.size(), array_Vertex.data());

		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
		glEnableVertexAttribArray(0);

		shaderYarnTess.SetParam(SHADER_PARAM_view_pos, camera.GetCameraPosition().x(), camera.GetCameraPosition().y(), camera.GetCameraPosition().z());
		shaderYarnTess.SetParam(SHADER_PARAM_view_dir, vd.x(), vd.y(), vd.z());
		shaderYarnTess.SetParam(SHADER_PARAM_light_pos_world, g_light_pos.x(), g_light_pos.y(), g_light_pos.z());
		ks::vec3 light_dir = g_light_pos.normalized();
		shaderYarnTess.SetParam(SHADER_PARAM_light_dir, light_dir.x(), light_dir.y(), light_dir.z());
		//before put the matrix into glsl, we have to transpose it because the diffenrence between row-major and column-major
		ks::mat4 gl_view = camera.GetMatrix();
		shaderYarnTess.SetParamMatrix4(SHADER_PARAM_view_matrix, gl_view.data());

		ks::mat4 gl_shadowMatrix = shadowMatrix;
		shaderYarnTess.SetParamMatrix4(SHADER_PARAM_shadow_matrix, gl_shadowMatrix.data());

		shaderYarnTess.SetParamMatrix4(SHADER_PARAM_shadow_R_matrix, shadowMat.data());

		shaderYarnTess.SetParam(SHADER_PARAM_Delta_0, mDelta_0);
		shaderYarnTess.SetParam(SHADER_PARAM_Delta_1, mDelta_1);
		shaderYarnTess.SetParam(SHADER_PARAM_Delta_2, mDelta_2);

		ks::mat4 gl_cameraMatrix_m4 = camera_matrix_m4;
		shaderYarnTess.SetParamMatrix4(SHADER_PARAM_camera_matrix, camera_matrix_m4.data());

		shaderYarnTess.SetParam(SHADER_PARAM_use_lod, fiberData.g_use_lod);
		shaderYarnTess.SetParam(SHADER_PARAM_use_lod_vis, fiberData.g_use_lod_vis);
		shaderYarnTess.SetParam(SHADER_PARAM_use_ao, fiberData.g_use_ao);
		shaderYarnTess.SetParam(SHADER_PARAM_use_diffuse, fiberData.g_use_diffuse);
		shaderYarnTess.SetParam(SHADER_PARAM_use_specular, fiberData.g_use_specular);
		shaderYarnTess.SetParam(SHADER_PARAM_use_regular_fiber, fiberData.g_use_regular_fiber);
		shaderYarnTess.SetParam(SHADER_PARAM_use_self_shadow, fiberData.g_use_self_shadow);
		shaderYarnTess.SetParam(SHADER_PARAM_use_shadow, fiberData.g_use_shadow);
		shaderYarnTess.SetParam(SHADER_PARAM_inv_ply_num, 1.0f / fiberData.g_ply_num);

		// specular
		shaderYarnTess.SetParam(SHADER_PARAM_color, fiberData.g_color[0], fiberData.g_color[1], fiberData.g_color[2]);
		shaderYarnTess.SetParam(SHADER_PARAM_lamda_R, fiberData.g_lamda_R);
		shaderYarnTess.SetParam(SHADER_PARAM_k_R, fiberData.g_k_R);

		// fiber parameters
		//std::cout << g_fiber_thickness << std::endl;
		shaderYarnTess.SetParam(SHADER_PARAM_fiber_thickness, fiberData.g_fiber_thickness);
		shaderYarnTess.SetParam(SHADER_PARAM_fiber_rho_min, fiberData.g_rho_min);
		shaderYarnTess.SetParam(SHADER_PARAM_fiber_rho_max, fiberData.g_rho_max);
		shaderYarnTess.SetParam(SHADER_PARAM_fiber_s_i, fiberData.g_s_i);
		shaderYarnTess.SetParam(SHADER_PARAM_inv_fiber_alpha, -1.0f / fiberData.g_alpha);
		//shaderYarnTess.SetParam(SHADER_PARAM_fiber_use_migration, g_use_migration);
		shaderYarnTess.SetParam(SHADER_PARAM_fiber_use_fly_away, fiberData.g_use_flyaways);

		// ply parameters
		shaderYarnTess.SetParam(SHADER_PARAM_inv_ply_alpha, -1.0f / fiberData.g_yarn_alpha);
		shaderYarnTess.SetParam(SHADER_PARAM_ply_radius, (fiberData.g_yarn_radius));
		shaderYarnTess.SetParam(SHADER_PARAM_ply_ellipse_long, fiberData.g_ellipse_long);
		shaderYarnTess.SetParam(SHADER_PARAM_ply_ellipse_short, fiberData.g_ellipse_short);
		shaderYarnTess.SetParam(SHADER_PARAM_ply_use_core_fibers, fiberData.g_use_core_fibers);

		shaderYarnTess.SetParam(SHADER_PARAM_scale, fiberData.g_scale);
		shaderYarnTess.SetParam(SHADER_PARAM_core_texture_height, float(g_core_texture_height));
		shaderYarnTess.SetParam(SHADER_PARAM_AD, AD.x(), AD.y(), AD.z());
		shaderYarnTess.SetParam(SHADER_PARAM_ATT, ATT.x(), ATT.y(), ATT.z());
		shaderYarnTess.SetParam(SHADER_PARAM_AR, AR.x(), AR.y(), AR.z());
		shaderYarnTess.SetParam(SHADER_PARAM_Beta_R, betaR);
		shaderYarnTess.SetParam(SHADER_PARAM_Beta_TT, betaTT);
		shaderYarnTess.SetParam(SHADER_PARAM_Beta_N, betaN);
		shaderYarnTess.SetParam(SHADER_PARAM_D_B, d_b);
		shaderYarnTess.SetParam(SHADER_PARAM_D_F, d_f);
		shaderYarnTess.SetParam(SHADER_PARAM_D_F_Inner, d_f_inner);
		shaderYarnTess.SetParam(SHADER_PARAM_D_Cross_Inter, d_cross_inter);
		shaderYarnTess.SetParam(SHADER_PARAM_Beta_Diffuse, beta_diffuse);

		// fly away
		shaderYarnTess.SetParam(SHADER_PARAM_fly_away_r_loop_max, fiberData.g_flyaway_loop_r1_mu / fiberData.g_ellipse_long);

		if (fiberData.g_fiber_num < 40.0)
			shaderYarnTess.SetParam(SHADER_PARAM_fly_away_rho_loop, abs(fiberData.g_flyaway_loop_density * fiberData.g_alpha / fiberData.g_fiber_num * 0.5f));
		else
			shaderYarnTess.SetParam(SHADER_PARAM_fly_away_rho_loop, abs(fiberData.g_flyaway_loop_density * fiberData.g_alpha / fiberData.g_fiber_num * 2.0f));

		// note 2.0 is a hard code number
		// it is difficult to determine the ratio of hair fibers
		if (fiberData.g_fiber_num < 40.0)
			shaderYarnTess.SetParam(SHADER_PARAM_fly_away_rho_hair, abs(fiberData.g_flyaway_hair_density * fiberData.g_alpha / fiberData.g_fiber_num * 1.0f));
		else
			shaderYarnTess.SetParam(SHADER_PARAM_fly_away_rho_hair, abs(fiberData.g_flyaway_hair_density * fiberData.g_alpha / fiberData.g_fiber_num * 2.0f));

		float g_flyaway_hair_rotation_scale = 1.0f / ((fiberData.g_flyaway_hair_ze_mu + fiberData.g_flyaway_hair_ze_sigma) * 0.5f) * fiberData.g_hair_rotation_scale;

		float g_flyaway_hair_length_scale = 1.0f / ((fiberData.g_flyaway_hair_pe_mu + fiberData.g_flyaway_hair_pe_sigma) * 0.25f) * fiberData.g_hair_length_scale
			* fiberData.g_yarn_radius / sqrt(fiberData.g_ellipse_long * fiberData.g_ellipse_short)
			* 0.03f * g_flyaway_hair_rotation_scale * 0.03f * g_flyaway_hair_rotation_scale * 0.03f * g_flyaway_hair_rotation_scale;

		//g_flyaway_hair_length_scale = g_flyaway_hair_length_scale * g_flyaway_hair_length_scale;

		//std::cout << g_flyaway_hair_rotation_scale << " " << g_flyaway_hair_length_scale << std::endl;

		shaderYarnTess.SetParam(SHADER_PARAM_fly_away_hair_rot_scale, 1.0f / g_flyaway_hair_rotation_scale);
		shaderYarnTess.SetParam(SHADER_PARAM_fly_away_hair_len_scale, g_flyaway_hair_length_scale);
		shaderYarnTess.SetParam(SHADER_PARAM_hair_scale, fiberData.g_hair_scale);

		shaderYarnTess.SetParam(SHADER_PARAM_ply_use_core_texture, fiberData.g_use_core_texture);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CROSS_SECTION_TEXUTRE);
		glBindTexture(GL_TEXTURE_1D, crossSectionTex);
		glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO_TEXTURE);
		glBindTexture(GL_TEXTURE_3D, CoreAOTex);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_LOD_BIAS, biasValue);

		//-----------------------	AO AZIMUTHAL TEXTURE -----------------------
		//8 directional AO, each as a 2D texture
		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO0_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[0]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO1_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[1]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO2_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[2]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO3_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[3]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO4_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[4]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO5_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[5]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO6_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[6]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO7_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs[7]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		//-----------------------	AO LONGITUDINAL TEXTURE -----------------------
		//8 directional AO, each as a 2D texture

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO0_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[0]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO1_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[1]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO2_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[2]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO3_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[3]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO4_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[4]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO5_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[5]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO6_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[6]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO7_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[7]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_AO8_TEXTURE_Long);
		glBindTexture(GL_TEXTURE_2D, CoreAOTexs_Long[8]);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);
		//--- end --------------

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SHADOW_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, depth_tex);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SHADOW_R_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, depth_debug_tex);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, coreTex);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_CORE_DIR_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, fbrDirTex);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glDrawArrays(GL_PATCHES, 0, vboYarnVertCount);

		if (fiberData.g_use_flyaways) {
			shaderFlyaway.BindProgram();

			shaderFlyaway.SetParamMatrix4(SHADER_PARAM_view_matrix, gl_view.data());
			shaderFlyaway.SetParamMatrix4(SHADER_PARAM_shadow_matrix, gl_shadowMatrix.data());
			shaderFlyaway.SetParamMatrix4(SHADER_PARAM_camera_matrix, camera_matrix_m4.data());
			shaderFlyaway.SetParam(SHADER_PARAM_light_dir, light_dir.x(), light_dir.y(), light_dir.z());
			shaderFlyaway.SetParam(SHADER_PARAM_view_dir, camera.GetCameraPosition().x(), camera.GetCameraPosition().y(), camera.GetCameraPosition().z());
			shaderFlyaway.SetParam(SHADER_PARAM_AD, AD.x(), AD.y(), AD.z());
			shaderFlyaway.SetParam(SHADER_PARAM_ATT, ATT.x(), ATT.y(), ATT.z());
			shaderFlyaway.SetParam(SHADER_PARAM_AR, AR.x(), AR.y(), AR.z());
			shaderFlyaway.SetParam(SHADER_PARAM_Beta_R, betaR);
			shaderFlyaway.SetParam(SHADER_PARAM_Beta_TT, betaTT);
			shaderFlyaway.SetParam(SHADER_PARAM_Beta_N, betaN);
			shaderFlyaway.SetParam(SHADER_PARAM_D_B, d_b);
			shaderFlyaway.SetParam(SHADER_PARAM_D_F, d_f);
			shaderFlyaway.SetParam(SHADER_PARAM_D_F_Inner, d_f_inner);
			shaderFlyaway.SetParam(SHADER_PARAM_D_Cross_Inter, d_cross_inter);
			shaderFlyaway.SetParam(SHADER_PARAM_Beta_Diffuse, beta_diffuse);
			shaderFlyaway.SetParam(SHADER_PARAM_TUBE_WIDTH, fiberData.g_yarn_radius);
			shaderFlyaway.SetParam(SHADER_PARAM_Delta_2, float(g_flyaway_texture_height - 10.0) / float(498));
			shaderFlyaway.SetParam(SHADER_PARAM_TEXTURE_OFFSET, mDelta_2);
			shaderFlyaway.SetParam(SHADER_PARAM_Delta_0, mDelta_0);
			shaderFlyaway.SetParam(SHADER_PARAM_flyaway_texture_height, float(g_flyaway_texture_height));

			glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SHADOW_TEXTURE);
			glBindTexture(GL_TEXTURE_2D, depth_tex);

			glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_FLYAWAY_TEXTURE);
			glBindTexture(GL_TEXTURE_2D, flyawayTex);
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

			glDrawArrays(GL_PATCHES, 0, vboYarnVertCount);
			glDisableVertexAttribArray(0);

		}

		/*std::string filename = "D:/test/imgui.png";
		saveTexture2DAsImage(SSAARenderBuffer->BufferTexID(), filename.data());*/

		//-------------------------------------------------------------- Object stage ---------------------------------------------------------

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, object_ibo);
		glBindBuffer(GL_ARRAY_BUFFER, object_vao);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
		glEnableVertexAttribArray(0);

		// specify normal attribute
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE,
			6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
		glEnableVertexAttribArray(1);

		shaderobject.BindProgram();

		shaderobject.SetParam(SHADER_PARAM_light_dir, light_dir.x(), light_dir.y(), light_dir.z());
		shaderobject.SetParamMatrix4(SHADER_PARAM_view_matrix, gl_view.data());
		shaderobject.SetParamMatrix4(SHADER_PARAM_shadow_matrix, gl_shadowMatrix.data());
		shaderobject.SetParam(SHADER_PARAM_view_pos, camera.GetCameraPosition().x(), camera.GetCameraPosition().y(), camera.GetCameraPosition().z());

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SHADOW_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, depth_tex);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, biasValue);

		glDrawElements(GL_TRIANGLES, object_tri_num * 3, GL_UNSIGNED_INT, NULL);

		SSAARenderBuffer->Unbind();

		//-------------------------------------------------------------- FINAL stage-----------------------------------------------------------
		//downsample

		glViewport(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		float vertices[] = {
			-1.0f, -1.0f, 0.5f,
			 1.0f, -1.0f, 0.5f,
			-1.0f,  1.0f, 0.5f,

			-1.0f,  1.0f, 0.5f,
			 1.0f, -1.0f, 0.5f,
			 1.0f,  1.0f, 0.5f,
		};

		unsigned int VAO, VBO;

		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);

		glBindVertexArray(VAO);

		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);

		shaderDownSample.BindProgram();
		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SSAA_TEXTURE);
		glBindTexture(GL_TEXTURE_2D, SSAARenderBuffer->BufferTexID());

		glActiveTexture(GL_TEXTURE0 + TEXTURE_UNIT_SSAA_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, SSAARenderBuffer2->BufferTexID());

		glDrawArrays(GL_TRIANGLES, 0, 6);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		//--------------------------------------------------------------  Start the Dear ImGui frame -----------------------------------------------------------

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		//----------------------------------------------------        light rotation/translation        ---------------------------------------------------------
		ImGui::Begin("Light Rotation");
		ImGui::SetWindowFontScale(2.f);
		ImGui::SetWindowSize(ImVec2(300, 650));
		lightControlDir = { camera.GetCameraPosition().x(), camera.GetCameraPosition().y(), camera.GetCameraPosition().z() };
		ImGui::gizmo3D("##modelGizmo", qRot, lightControlDir, 250, imguiGizmo::modeDual | imguiGizmo::cubeAtOrigin);
		
		vgm::Mat4 modelMat4 = mat4_cast(qRot);
		vgm::Vec4 lightdir = modelMat4 * vgm::Vec4(0, 0, 1, 1.f);
		g_light_pos = ks::vec3(lightdir.x, lightdir.y, lightdir.z) * 25;

		if (ImGui::Button("Reset"))
		{
			qRot = vgm::Quat(1.f, 0.f, 0.f, 0.f);
		}

		//----------------------------------------------------        camera rotation/translation        ---------------------------------------------------------
		ImGuiIO& io = ImGui::GetIO();

		float wheel_delta = io.MouseWheel;
		if (wheel_delta != 0.0f) {
			camera.AddDistance(-wheel_delta * 0.4);
		}

		ImVec2 mouse_pos = ImGui::GetMousePos();
		ImVec2 diff = ImVec2(mouse_pos.x - last_mouse_pos.x, mouse_pos.y - last_mouse_pos.y);
		if ((mouse_pos.x != last_mouse_pos.x) || (mouse_pos.y != last_mouse_pos.y)) {
			last_mouse_pos = mouse_pos;
		}

		if (ImGui::IsMouseDown(1)) 
		{
				camera.AddRotation(diff.x * 0.01f, diff.y * 0.01f);
		}
		 

		if (ImGui::SliderFloat("yarn_thickness", &thickness, 0.0f, 1.f)) {
			fiberData.g_yarn_radius = thickness * 0.1;
		}

		ImGui::SliderFloat("LOD_Bias", &biasValue, -5, 5);
		
		bool flyaway_regenerate = false;
		if (ImGui::SliderFloat("hair_density", &fiberData.g_flyaway_hair_density, 0, 120)) {
			flyaway_regenerate = true;
		}
		if (ImGui::SliderFloat("loop_density", &fiberData.g_flyaway_loop_density, 0, 120)) {
			flyaway_regenerate = true;
		}
		if (ImGui::SliderFloat("loop_radius", &fiberData.g_flyaway_loop_r1_mu, 0, 0.1)) {
			flyaway_regenerate = true;
		}
		if (ImGui::SliderFloat("hair_radius", &fiberData.g_flyaway_hair_re_mu, 0, 0.1)) {
			flyaway_regenerate = true;
		}
		if (flyaway_regenerate) UpdateOfflineYarn();

		ImGui::End();
		ImGui::Render();
		// end for rotation control panel
		//double currentTime = glfwGetTime();
		//nbFrames++;
		//if (currentTime - lastTime >= 1.0) { // If last prinf() was more than 1 sec ago
		//	// printf and reset timer
		//	printf("%f ms/frame\n", 1000.0 / double(nbFrames));
		//	nbFrames = 0;
		//	lastTime += 1.0;
		//}

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
}