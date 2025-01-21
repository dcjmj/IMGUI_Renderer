#pragma once
#include <windows.h>
#include <GL/glew.h>      // GLEW must be included before any OpenGL header
#include <GLFW/glfw3.h>   // GLFW includes OpenGL headers internally
#include <GL/GLU.h>       // Optional, depends on your usage
#include <glm/glm.hpp>    // GLM is a math library, safe to include last
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 

#include "imGuIZMOquat.h"
#include "imgui.h"
#include "../imgui/backends/imgui_impl_glfw.h"
#include "../imgui/backends/imgui_impl_opengl3.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "maths.h"	
#include "fiber.h"
#include "ShaderTypes.h"
#include "camera.h"
#include "globals.h"

#define			ATTRIB_VERTEX			0
#define			ATTRIB_TEXTURE			1
#define			ATTRIB_NORMAL			2
#define			ATTRIB_COLOR			3
#define			DEFAULT_WIDTH			1920//1317//1317/1280//1920//1280
#define			DEFAULT_HEIGHT			1280//1317/960/1080//960
#define			CHECKER_BOARD_WIDTH		1024
#define			CHECKER_BOARD_HEIGHT	1024
#define			SHADOW_MAP_SIZE			4096

#define CY_GLSL_INVALID_ID 0xFFFFFFFF

static const GLfloat	ones[] = { 1.0f };
static const GLfloat	zero[] = { 0.0f };
static const GLfloat	gray[] = { 1.0f, 1.0f, 1.0f, 1.0f };

#define PI	3.1415926535897932384626433832795f

#define SWAP(type,a,b) { type _t=(a); (a)=(b); (b)=_t; }

#define a2rad(a) ((a)*PI/180.0f)
class OpenGLRenderBuffer
{
public:
	enum Format {
		FORMAT_INT_8,
		FORMAT_INT_16,
		FORMAT_INT_32,
		FORMAT_UINT_8,
		FORMAT_UINT_16,
		FORMAT_UINT_32,
		FORMAT_FLOAT_16,
		FORMAT_FLOAT_32,
	};

	OpenGLRenderBuffer() : buffer(0), bufferDepth(0) {}

	void Initialize(bool useDepthBuffer, int num_colorbuffer = 1);
	void Resize(int width, int height, Format format, int numChannels);
	void Delete();

	void Bind();
	void Unbind();

	unsigned int BufferTexID() { return bufferTex; }

	unsigned int buffer, bufferTex, bufferDepth, buffer2Tex;
	int prevBuffer;
};

//GLSL SHADER Class

/// Shader class
class cyGLSLShader
{
private:
	GLuint shaderID;
public:
	cyGLSLShader() : shaderID(CY_GLSL_INVALID_ID) {}
	virtual ~cyGLSLShader() { Clear(); }

	GLuint ShaderID() const { return shaderID; }

	void Clear()
	{
		if (shaderID != CY_GLSL_INVALID_ID) {
			glDeleteShader(shaderID);
			shaderID = CY_GLSL_INVALID_ID;
		}
	}

	bool CompileShaderFile(const char* filename, GLenum shaderType, std::ostream& outStream = std::cout)
	{
		std::ifstream shaderStream(filename, std::ios::in);
		if (!shaderStream.is_open()) {
			outStream << "ERROR: Cannot open file." << std::endl;
			return false;
		}

		std::string shaderSourceCode((std::istreambuf_iterator<char>(shaderStream)), std::istreambuf_iterator<char>());
		shaderStream.close();

		return CompileShaderSource(shaderSourceCode.data(), shaderType, outStream);
	}

	bool CompileShaderSource(const char* shaderSourceCode, GLenum shaderType, std::ostream& outStream = std::cout)
	{
		Clear();

		shaderID = glCreateShader(shaderType);
		glShaderSource(shaderID, 1, &shaderSourceCode, NULL);
		glCompileShader(shaderID);

		GLint result = GL_FALSE;
		glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);

		int infoLogLength;
		glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
		if (infoLogLength > 1) {
			std::vector<char> compilerMessage(infoLogLength);
			glGetShaderInfoLog(shaderID, infoLogLength, NULL, compilerMessage.data());
			outStream << "ERROR: " << compilerMessage.data() << std::endl;
		}

		if (result) {
			GLint stype;
			glGetShaderiv(shaderID, GL_SHADER_TYPE, &stype);
			if (stype != shaderType) {
				outStream << "ERROR: Incorrect shader type." << std::endl;
				return false;
			}
		}

		return result == GL_TRUE;
	}
};

class cyGLSLProgram
{
private:
	GLuint programID;
	std::vector<GLint> params;
public:
	cyGLSLProgram() : programID(CY_GLSL_INVALID_ID) {}
	virtual ~cyGLSLProgram() { Clear(); }

	GLuint ProgramID() const { return programID; }

	void Clear()
	{
		if (programID != CY_GLSL_INVALID_ID) {
			glDeleteProgram(programID);
			programID = CY_GLSL_INVALID_ID;
		}
	}

	void CreateProgram() { Clear(); programID = glCreateProgram(); }

	void BindProgram() { glUseProgram(programID); }

	void AttachShader(const cyGLSLShader& shader) { AttachShader(shader.ShaderID()); }
	void AttachShader(GLuint shaderID) { glAttachShader(programID, shaderID); }

	bool LinkProgram(std::ostream& outStream = std::cout)
	{
		glLinkProgram(programID);

		GLint result = GL_FALSE;
		glGetProgramiv(programID, GL_LINK_STATUS, &result);

		int infoLogLength;
		glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
		if (infoLogLength > 1) {
			std::vector<char> compilerMessage(infoLogLength);
			glGetProgramInfoLog(programID, infoLogLength, NULL, compilerMessage.data());
			outStream << "ERROR: " << compilerMessage.data() << std::endl;
		}

		return result == GL_TRUE;
	}

	bool BuildProgramFiles(const char* vertexShaderFile,
		const char* fragmentShaderFile = NULL,
		const char* geometryShaderFile = NULL,
		const char* tessControlShaderFile = NULL,
		const char* tessEvaluationShaderFile = NULL,
		std::ostream& outStream = std::cout)
	{
		Clear();
		CreateProgram();
		cyGLSLShader vs, fs, gs, tcs, tes;
		std::stringstream output;
		if (!vs.CompileShaderFile(vertexShaderFile, GL_VERTEX_SHADER, output)) {
			outStream << "ERROR: Failed compiling vertex shader \"" << vertexShaderFile << ".\"" << std::endl << output.str();
			return false;
		}
		AttachShader(vs);
		if (fragmentShaderFile) {
			if (!fs.CompileShaderFile(fragmentShaderFile, GL_FRAGMENT_SHADER, output)) {
				outStream << "ERROR: Failed compiling fragment shader \"" << fragmentShaderFile << ".\"" << std::endl << output.str();
				return false;
			}
			AttachShader(fs);
		}
		if (geometryShaderFile) {
			if (!gs.CompileShaderFile(geometryShaderFile, GL_GEOMETRY_SHADER, output)) {
				outStream << "ERROR: Failed compiling geometry shader \"" << geometryShaderFile << ".\"" << std::endl << output.str();
				return false;
			}
			AttachShader(gs);
		}
		if (tessControlShaderFile) {
			if (!tcs.CompileShaderFile(tessControlShaderFile, GL_TESS_CONTROL_SHADER, output)) {
				outStream << "ERROR: Failed compiling tessellation control shader \"" << tessControlShaderFile << ".\"" << std::endl << output.str();
				return false;
			}
			AttachShader(tcs);
		}
		if (tessEvaluationShaderFile) {
			if (!tes.CompileShaderFile(tessEvaluationShaderFile, GL_TESS_EVALUATION_SHADER, output)) {
				outStream << "ERROR: Failed compiling tessellation evaluation shader \"" << tessEvaluationShaderFile << ".\"" << std::endl << output.str();
				return false;
			}
			AttachShader(tes);
		}
		LinkProgram(outStream);
		return true;
	}

	bool BuildProgramSources(const char* vertexShaderSourceCode,
		const char* fragmentShaderSourceCode,
		const char* geometryShaderSourceCode = NULL,
		const char* tessControlShaderSourceCode = NULL,
		const char* tessEvaluationShaderSourceCode = NULL,
		std::ostream& outStream = std::cout)
	{
		Clear();
		CreateProgram();
		cyGLSLShader vs, fs, gs, tcs, tes;
		std::stringstream output;
		if (!vs.CompileShaderSource(vertexShaderSourceCode, GL_VERTEX_SHADER, output)) {
			outStream << "ERROR: Failed compiling vertex shader." << std::endl << output.str();
			return false;
		}
		AttachShader(vs);
		if (!fs.CompileShaderSource(fragmentShaderSourceCode, GL_FRAGMENT_SHADER, output)) {
			outStream << "ERROR: Failed compiling fragment shader." << std::endl << output.str();
			return false;
		}
		AttachShader(fs);
		if (geometryShaderSourceCode) {
			if (!gs.CompileShaderSource(geometryShaderSourceCode, GL_GEOMETRY_SHADER, output)) {
				outStream << "ERROR: Failed compiling geometry shader." << std::endl << output.str();
				return false;
			}
			AttachShader(gs);
		}
		if (tessControlShaderSourceCode) {
			if (!tcs.CompileShaderSource(tessControlShaderSourceCode, GL_TESS_CONTROL_SHADER, output)) {
				outStream << "ERROR: Failed compiling tessellation control shader." << std::endl << output.str();
				return false;
			}
			AttachShader(tcs);
		}
		if (tessEvaluationShaderSourceCode) {
			if (!tes.CompileShaderSource(tessEvaluationShaderSourceCode, GL_TESS_EVALUATION_SHADER, output)) {
				outStream << "ERROR: Failed compiling tessellation evaluation shader." << std::endl << output.str();
				return false;
			}
			AttachShader(tes);
		}
		LinkProgram(outStream);
		return true;
	}

	/// Registers a single parameter.
	/// The id must be unique and the name should match a uniform parameter name in one of the shaders.
	/// The id values for different parameters don't have to be consecutive, but unused id values waste memory.
	void RegisterParam(unsigned int id, const char* name, std::ostream& outStream = std::cout)
	{
		if (params.size() <= id) params.resize(id + 1, -1);
		params[id] = glGetUniformLocation(programID, name);
		if (params[id] < 0) {
			GLenum error = glGetError();
			GLenum newError;
			while ((newError = glGetError()) != GL_NO_ERROR) error = newError; // get the latest error.
			outStream << "OpenGL ERROR: " << gluErrorString(error) << ". Parameter \"" << name << "\" could not be registered." << std::endl;
		}
	}

	/// Registers multiple parameters.
	/// The names should be separated by a space character.
	void RegisterParams(const char* names, unsigned int startingID = 0, std::ostream& outStream = std::cout)
	{
		std::stringstream ss(names);
		unsigned int id = startingID;
		while (ss.good()) {
			std::string name;
			ss >> name;
			RegisterParam(id++, name.c_str(), outStream);
		}
	}

	void SetParam(unsigned int paramID, float x) { glUniform1f(params[paramID], x); }
	void SetParam(unsigned int paramID, float x, float y) { glUniform2f(params[paramID], x, y); }
	void SetParam(unsigned int paramID, float x, float y, float z) { glUniform3f(params[paramID], x, y, z); }
	void SetParam(unsigned int paramID, float x, float y, float z, float w) { glUniform4f(params[paramID], x, y, z, w); }
	void SetParam(unsigned int paramID, int x) { glUniform1i(params[paramID], x); }
	void SetParam(unsigned int paramID, int x, int y) { glUniform2i(params[paramID], x, y); }
	void SetParam(unsigned int paramID, int x, int y, int z) { glUniform3i(params[paramID], x, y, z); }
	void SetParam(unsigned int paramID, int x, int y, int z, int w) { glUniform4i(params[paramID], x, y, z, w); }
	void SetParam(unsigned int paramID, unsigned int x) { glUniform1ui(params[paramID], x); }
	void SetParam(unsigned int paramID, unsigned int x, unsigned int y) { glUniform2ui(params[paramID], x, y); }
	void SetParam(unsigned int paramID, unsigned int x, unsigned int y, unsigned int z) { glUniform3ui(params[paramID], x, y, z); }
	void SetParam(unsigned int paramID, unsigned int x, unsigned int y, unsigned int z, unsigned int w) { glUniform4ui(params[paramID], x, y, z, w); }
	void SetParams(unsigned int paramID, unsigned int count, const float* m) { glUniform1fv(params[paramID], count, m); }
	void SetParamsV3(unsigned int paramID, unsigned int count, const float* m) { glUniform3fv(params[paramID], count, m); }
	void SetParamMatrix2(unsigned int paramID, const float* m) { glUniformMatrix2fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix2x3(unsigned int paramID, const float* m) { glUniformMatrix2x3fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix2x4(unsigned int paramID, const float* m) { glUniformMatrix2x4fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix3x2(unsigned int paramID, const float* m) { glUniformMatrix3x2fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix3(unsigned int paramID, const float* m) { glUniformMatrix3fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix3x4(unsigned int paramID, const float* m) { glUniformMatrix3x4fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix4x2(unsigned int paramID, const float* m) { glUniformMatrix4x2fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix4x3(unsigned int paramID, const float* m) { glUniformMatrix4x3fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamMatrix4(unsigned int paramID, const float* m) { glUniformMatrix4fv(params[paramID], 1, GL_FALSE, m); }
	void SetParamsMatrix4(unsigned int paramID, unsigned int count, const float* m) { glUniformMatrix4fv(params[paramID], 16, GL_FALSE, m); }

};

struct Face {
	int vertexIndices[3];
	int normalIndices[3];
};

struct ObjData {
	std::vector<ks::vec3> vertices;
	std::vector<ks::vec2> texCoords;
	std::vector<ks::vec3> normals;
	std::vector<Face> faces;
};

class GLWidget {
public:
	void render_task(float* Yarn_ctrPoints, int* first_ctrP_idx, int yarn_num, Globals &globals);
	void update_yarn_buffer(float* Yarn_ctrPoints, int* first_ctrP_idx, int yarn_num);
	void initializeGL();
	void UpdateOfflineYarn();
	void UpdateCrossSectionTexture();
	void LoadCoreDirImage();
	void LoadCoreImage();
	void InitOfflineYarns();
	void Init1DCylinderNormTexture();
	void InitCoreAOTexture();
	void InitCoreAOTexture_Long();
	bool InitShaders();
	void InitObject(std::string filepath = "");
	void InitObject(const Mesh &mesh); 
	void UpdateShadowTexture();
	void UpdateConfig2OfflineYarn(Fiber::Yarn* y);
	void ReadFrameBuffer();
	void CreateIndexBuffer(const ObjData& objData);
	void CreateIndexBuffer(const Mesh &mesh);
	GLuint defaultFramebufferObject() const;
	void LoadIntegration(char* filename);
	void UpdateEnvShadowTexture();
	void InitEnvShadowTextureBuffer();
	void InitShadowTextureBuffer();
	void createRandomTexture(int width, int height);

	GLuint object_tri_num;

	GLFWwindow* window;

	GLuint object_vao;
	GLuint object_ibo;

	GLuint cloth_vao;
	GLuint cloth_ibo;
	GLuint RandomTex;

	std::vector<float> lambdas_final; // Initial lambda values
	std::vector<ks::vec3> centers_final;
	Eigen::VectorXf weights_final;
	Eigen::VectorXf weights;
	bool g_use_envmap = false;


	GLuint WINDOW_WIDTH = DEFAULT_WIDTH;
	GLuint WINDOW_HEIGHT = DEFAULT_HEIGHT;


	int                                     g_integral_width;
	int                                     g_integral_height;
	int                                     g_integral_angle_res;
	int                                     g_integral_lambda_res;

	//env map shadow
	GLuint          depth_fbo_env[16];
	GLuint			depth_tex_env[16];
	ks::mat4								light_proj_matrix_env[16];
	ks::mat4								light_view_matrix_env[16];
	ks::mat4								light_matrix_env[16];
	ks::vec3 envlight_pos[16];

	ks::vec4								lightPosCam4;
	ks::mat4							    lightProjMatrix;
	ks::mat4								lightViewMatrix;

	ks::mat4								mIdentity;
	ks::mat4								mI4;
	ks::mat4								light_proj_matrix;
	ks::mat4								light_view_matrix;
	ks::mat4								camera_matrix_m4;
	ks::vec3								adjust_color;
	ks::vec3								residual_color[8];

	ks::Camera								camera;
	ks::vec3 center;

	//object shadow
	GLuint          depth_fbo;
	//	QOpenGLFramebufferObject* depth_fbo;
	GLuint          depth_tex;
	GLuint          depth_debug_tex;
	GLuint									IntegralTex;
	std::vector<GLfloat> array_Vertex;
	bool fiber_updated = true;

	//related setup
	int g_model_idx = 1;
	int g_config_idx = 1;
	int g_shading_idx = 1;

	//self AO texture
	ks::vec3 AR = ks::vec3(0.15, 0.15, 0.15);
	ks::vec3 AD = ks::vec3(0.39, 0.305, 0.04);
	ks::vec3 ATT = ks::vec3(0.3, 0.21, 0.11);
	float betaR = 0.2;
	float betaTT = 0.2;
	float betaN = 0.6;
	float d_b = 0.8;
	float d_f = 0.8;
	float beta_diffuse = 0.90;
	float d_f_inner = 0.8;
	float d_cross_inter = 0.0;

	FiberGenerationData fiberData;
	FiberGenerationData* Get_FiberData() { return &fiberData; }

	Fiber::Yarn* offlinePlys;
	Fiber::Yarn* offlineFlyaway;
	Fiber::Yarn* offlineYarns;
	Fiber::Yarn* offlineCore;
	std::vector<std::vector<ks::vec3>>		offlineYarnVertices;
	std::vector<std::vector<ks::vec3>>		offlineCoreVertices;
	GLuint* offlineYarnBuffer;
	std::vector<GLuint> offlineCoreBuffer;				// used for core texture generation
	OpenGLRenderBuffer* offlineCoreRenderBuffer;		// write core texture to this buffer
	OpenGLRenderBuffer* offlineCoreTangentRenderBuffer;		// write core tangent texture to this buffer
	OpenGLRenderBuffer* SSAARenderBuffer;		// write core tangent texture to this buffer
	OpenGLRenderBuffer* SSAARenderBuffer2;		// write core tangent texture to this buffer

	int										g_core_texture_width;
	int										g_core_texture_height;

	int										g_flyaway_texture_width;
	int										g_flyaway_texture_height;

	int										g_core_AO_texture_width;
	int										g_core_AO_texture_height;

	int										g_flyaway_texture_left;
	int										g_flyaway_texture_right;
	int										g_flyaway_texture_top = DEFAULT_HEIGHT / 2;
	int										g_flyaway_texture_bot = DEFAULT_HEIGHT / 2;
	int										g_regular_width;

	GLint									DefaultShaderProgram;
	GLuint									ttIntegralTex;
	GLuint									rIntegralTex;
	GLuint									cylinderTex;
	GLuint									crossSectionTex;
	GLuint									checkerBoardTex;
	GLuint									coreTex;
	GLuint									flyawayTex;
	GLuint									fbrDirTex;
	GLuint									selfShadowTex;
	GLuint									indirectLightTex;
	GLuint									ttLightTex;
	GLuint									DensityTex;
	GLuint									AOTex;
	GLuint 								    CoreAOTex;
	GLuint 								    CoreAOTexs[8];
	GLuint 								    CoreAOTexs_Long[9];
	GLuint									vertexArrays[VERTEX_ARRAY_COUNT];
	GLuint									vboYarn[VBO_YARN_COUNT];
	float biasValue = -10;
	float mDelta_0 = 0.14;
	float mDelta_1 = 1;
	float mDelta_2 = 1;

	std::vector<int>						vboYarnVertStart;
	int										vboYarnVertCount;

	GLuint									objectTex;

	float					g_center_offset;
	//ks::vec3				g_light_pos = ks::vec3(0,0,-50);
	ks::vec3				g_light_pos = ks::vec3(10, 10, 4);
	float			fov = 60.0f;
	float			aspect = 1.0f;
	float			znear = 0.01f;
	float			zfar = 100.0f;

	//ALL SHADER PROGRAM
	cyGLSLProgram		shaderYarnTess, shaderCheckerBoard, shaderYarnVerts,
		shaderOfflineYarn, shaderOfflineCore, shaderOfflineFbrDir, shaderFrameBuffer,
		shaderDepthBuffer, shaderYarnTubeDepth, shaderYarnTube, shaderObject, shaderObjectDepth, shaderTest, shaderOfflineCoreT, shaderFlyaway, shaderDownSample, shaderobjectshadow, shaderobject;
};