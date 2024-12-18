#pragma once
#include"maths.h"


enum {
	SHADER_PARAM_model_matrix,
	SHADER_PARAM_proj_matrix,
	SHADER_PARAM_view_matrix,
	SHADER_PARAM_camera_matrix,
	SHADER_PARAM_shadow_matrix,
	SHADER_PARAM_shadow_R_matrix,

	SHADER_PARAM_view_dir,
	SHADER_PARAM_view_pos,

	SHADER_PARAM_light_pos,
	SHADER_PARAM_light_pos_world,
	SHADER_PARAM_light_dir,

	SHADER_PARAM_color,
	SHADER_PARAM_lamda_R,
	SHADER_PARAM_k_R,

	SHADER_PARAM_checker_mv_matrix,
	SHADER_PARAM_checker_proj_matrix,
	SHADER_PARAM_checker_texture,

	SHADER_PARAM_cross_section_texture,
	SHADER_PARAM_cylinder_texture,
	SHADER_PARAM_frame_buffer_texture,
	SHADER_PARAM_depth_buffer_texture,
	SHADER_PARAM_self_shadow_texture,
	SHADER_PARAM_shadow_texture,
	SHADER_PARAM_flyaway_texture,
	SHADER_PARAM_core_texture,
	SHADER_PARAM_core_dir_texture,
	SHADER_PARAM_core_AO_texture,
	SHADER_PARAM_core_AO0_texture,
	SHADER_PARAM_core_AO1_texture,
	SHADER_PARAM_core_AO2_texture,
	SHADER_PARAM_core_AO3_texture,
	SHADER_PARAM_core_AO4_texture,
	SHADER_PARAM_core_AO5_texture,
	SHADER_PARAM_core_AO6_texture,
	SHADER_PARAM_core_AO7_texture,
	SHADER_PARAM_core_AO0_texture_Long,
	SHADER_PARAM_core_AO1_texture_Long,
	SHADER_PARAM_core_AO2_texture_Long,
	SHADER_PARAM_core_AO3_texture_Long,
	SHADER_PARAM_core_AO4_texture_Long,
	SHADER_PARAM_core_AO5_texture_Long,
	SHADER_PARAM_core_AO6_texture_Long,
	SHADER_PARAM_core_AO7_texture_Long,
	SHADER_PARAM_core_AO8_texture_Long,
	SHADER_PARAM_SSAA_texture,
	SHADER_PARAM_SSAA_texture2,
	SHADER_PARAM_shadow_R_texture,
	SHADER_PARAM_TUBE_WIDTH,
	SHADER_PARAM_Light_Dir,

	SHADER_PARAM_fiber_thickness,
	SHADER_PARAM_fiber_rho_min,
	SHADER_PARAM_fiber_rho_max,
	SHADER_PARAM_fiber_s_i,
	SHADER_PARAM_inv_fiber_alpha,
	//SHADER_PARAM_fiber_use_migration,
	SHADER_PARAM_fiber_use_fly_away,
	SHADER_PARAM_ply_ellipse_long,
	SHADER_PARAM_ply_ellipse_short,
	SHADER_PARAM_inv_ply_alpha,
	SHADER_PARAM_ply_radius,
	SHADER_PARAM_inv_ply_num,
	SHADER_PARAM_fly_away_r_loop_max,
	SHADER_PARAM_fly_away_rho_loop,
	SHADER_PARAM_fly_away_rho_hair,
	SHADER_PARAM_fly_away_hair_rot_scale,
	SHADER_PARAM_fly_away_hair_len_scale,
	SHADER_PARAM_ply_use_core_fibers,
	SHADER_PARAM_ply_use_core_texture,

	SHADER_PARAM_scale,

	SHADER_PARAM_use_lod,
	SHADER_PARAM_use_ao,
	SHADER_PARAM_use_diffuse,
	SHADER_PARAM_use_specular,
	SHADER_PARAM_use_regular_fiber,
	SHADER_PARAM_use_shadow,
	SHADER_PARAM_use_self_shadow,
	SHADER_PARAM_use_lod_vis,
	SHADER_PARAM_center_offset,
	SHADER_PARAM_core_texture_height,
	SHADER_PARAM_flyaway_texture_height,
	SHADER_PARAM_hair_scale,

	SHADER_PARAM_AD,
	SHADER_PARAM_AR,
	SHADER_PARAM_ATT,
	SHADER_PARAM_Beta_R,
	SHADER_PARAM_Beta_TT,
	SHADER_PARAM_Beta_N,
	SHADER_PARAM_D_B,
	SHADER_PARAM_D_F,
	SHADER_PARAM_Beta_Diffuse,
	SHADER_PARAM_D_F_Inner,
	SHADER_PARAM_D_Cross_Inter,
	SHADER_PARAM_DENSITY_TEXTURE,
	SHADER_PARAM_GRID_RES,
	SHADER_PARAM_GRID_RADIUS,
	SHADER_PARAM_SCENE_MIN,
	SHADER_PARAM_Delta_0,
	SHADER_PARAM_Delta_1,
	SHADER_PARAM_Delta_2,
	SHADER_PARAM_TEXTURE_OFFSET,
};

//-------------------------------------------------------------------------------
// other enum

enum {
	TEXTURE_UNIT_CHECKER_TEXUTRE,
	TEXTURE_UNIT_CROSS_SECTION_TEXUTRE,
	TEXTURE_UNIT_CYLINDER_TEXUTRE,
	TEXTURE_UNIT_FRAMEBUFFER_TEXUTRE,
	TEXTURE_UNIT_DEPTHBUFFER_TEXUTRE,
	TEXTURE_UNIT_SELFSHADOW_TEXTURE,
	TEXTURE_UNIT_INDIRECT_LIGHT_TEXTURE,
	TEXTURE_UNIT_SHADOW_TEXTURE,
	TEXTURE_UNIT_FLYAWAY_TEXTURE,
	TEXTURE_UNIT_CORE_TEXTURE,
	TEXTURE_UNIT_AO_TEXTURE,
	TEXTURE_UNIT_FIBER_DIR_TEXTURE,
	TEXTURE_UNIT_TT_LIGHT_TEXTURE,
	TEXTURE_UNIT_OBJ_AMB_TEXTURE,
	TEXTURE_UNIT_R_INTEGRAL_TEXTURE,
	TEXTURE_UNIT_TT_INTEGRAL_TEXTURE,
	TEXTURE_UNIT_DENSITY_TEXTURE,
	TEXTURE_UNIT_CORE_DIR_TEXTURE,
	TEXTURE_UNIT_CORE_AO_TEXTURE,
	TEXTURE_UNIT_CORE_AO0_TEXTURE,
	TEXTURE_UNIT_CORE_AO1_TEXTURE,
	TEXTURE_UNIT_CORE_AO2_TEXTURE,
	TEXTURE_UNIT_CORE_AO3_TEXTURE,
	TEXTURE_UNIT_CORE_AO4_TEXTURE,
	TEXTURE_UNIT_CORE_AO5_TEXTURE,
	TEXTURE_UNIT_CORE_AO6_TEXTURE,
	TEXTURE_UNIT_CORE_AO7_TEXTURE,
	TEXTURE_UNIT_CORE_AO0_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO1_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO2_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO3_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO4_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO5_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO6_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO7_TEXTURE_Long,
	TEXTURE_UNIT_CORE_AO8_TEXTURE_Long,
	TEXTURE_UNIT_SSAA_TEXTURE,
	TEXTURE_UNIT_SSAA_TEXTURE2,
	TEXTURE_UNIT_SHADOW_R_TEXTURE,
};

enum {
	VERTEX_ARRAY_YARN_CONTROL_POINTS,
	VERTEX_ARRAY_COUNT
};

enum {
	VBO_YARN_VERTEX
};

# define VBO_YARN_COUNT (VBO_YARN_VERTEX+1)

typedef struct FiberGenerationData {
	ks::vec3				g_aabb_micro_ct_pMin = ks::vec3(-0.0498028, -0.0498028, -0.486193);
	float 					g_ellipse_long = 0.0282957;
	ks::vec3 				g_aabb_micro_ct_pMax = ks::vec3(0.0498028, 0.0498028, 0.486193);
	float 					g_ellipse_short = 0.0208841;
	ks::vec3				g_color = ks::vec3(1.0f, 0.8f, 0.1f);
	float 					g_yarn_radius = 0.02866;
	float 					g_fiber_thickness = 0.008f;
	float 					g_yarn_alpha = 0.38;
	float 					g_fiber_num = 75;
	float 					g_ply_num = 2;
	float					g_direct_light_intensity = 0.6f;
	float					g_amb_light_intensity = 0.5f;
	float					g_center_offset;
	float 					g_scale = 1.0f;
	float					g_hair_length_scale = 1.f;
	float					g_hair_rotation_scale = 1.0f;
	float					g_hair_scale = 1.7f;
	float					g_lamda_R = 0.2f;
	float					g_k_R = 0.5f;


	//-------------------------------------------------------------------------------
	// fly away parameters
	float 					g_flyaway_loop_density = 22.17;
	float 					g_flyaway_hair_density = 33.7667;
	float 					g_flyaway_hair_ze_mu = -0.00252956;
	float 					g_flyaway_hair_ze_sigma = 0.0572919;
	float 					g_flyaway_hair_r0_mu = 0.0197073;
	float 					g_flyaway_hair_r0_sigma = 0.00556197;
	float 					g_flyaway_hair_re_mu = 0.0162684;
	float 					g_flyaway_hair_re_sigma = 0.00850055;
	float 					g_flyaway_hair_pe_mu = 0.376521;
	float 					g_flyaway_hair_pe_sigma = 0.325502;
	float 					g_flyaway_loop_r1_mu = 0.0237306;
	float 					g_flyaway_loop_r1_sigma = 0.00537727;
	float 					g_z_step_size = 0.01;
	float 					g_z_step_num = 98;
	float                   g_epsilon = 0.0521;
	float                   g_r_max = 1;
	float                   g_beta = 0.292756;
	float                   g_alpha = 0.38;
	float                   g_s_i = 1.0;
	float                   g_rho_min = 0.1;
	float                   g_rho_max = 0.85;

	bool 					g_yarn_clock_wise = false;
	bool 					g_use_flyaways = true;
	bool 					g_use_migration = true;
	bool                    g_fiber_clock_wise = true;

	//Render strategy
	bool 					g_use_core_fibers = true;
	bool 					g_use_core_texture = true;
	bool 					g_use_ao = true;
	bool 					g_use_direct = true;
	bool 					g_use_iso = true;
	bool 					g_use_shadow = true;
	bool 					g_use_self_shadow = true;
	bool 					g_use_regular_fiber = true;
	bool 					g_use_lod = true;
	bool 					g_use_lod_vis = false;
	bool 					g_draw_light_source = false;
	bool					g_draw_reference = false;
	bool					g_use_phys_shading = true;
	bool					g_draw_checkerboard = false;
	bool					g_draw_object = false;
	bool					g_use_precomputed_amb = true;
	bool 					g_use_diffuse = true;
	bool 					g_use_specular = true;

} FiberGenerationData;