 /*
 Copyright (c) 2015 to 2016 by Cornell University and The Regents Of
 The University Of California. All Rights Reserved.

 Permission to use this Procedural Yarn Fitting and Generation Tool (the "Work")
 and its associated copyrights solely for educational, research and non-profit
 purposes, without fee is hereby granted, provided that the user agrees as
 follows:

 Those desiring to incorporate the Work into commercial products or use Work and
 its associated copyrights for commercial purposes should contact the Center for
 Technology Licensing at Cornell University at

 395 Pine Tree Road, Suite 310, Ithaca, NY 14850;
 email: ctl-connect@cornell.edu;
 Tel: 607-254-4698;
 FAX: 607-254-5454

 for a commercial license.

 IN NO EVENT SHALL CORNELL UNIVERSITY ("CORNELL") OR THE UNIVERSITY OF
 CALIFORNIA ("UC") BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
 INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
 THE USE OF THE WORK AND ITS ASSOCIATED COPYRIGHTS, EVEN IF CORNELL OR UC MAY
 HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 THE WORK PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND NEITHER CORNELL NOR UC HAS
 ANY OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 MODIFICATIONS. CORNELL AND UC MAKE NO REPRESENTATIONS AND EXTEND NO WARRANTIES
 OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR
 THAT THE USE OF WORK AND ITS ASSOCIATED COPYRIGHTS WILL NOT INFRINGE ANY PATENT,
 TRADEMARK OR OTHER RIGHTS.
 */

 #pragma once
 //#define VERBOSE
 #define IMPROVED_FLYAWAYS
 #include <omp.h>
 #include "maths.h"
 #include <iostream>
 #include <fstream>
 #include "aabb.h"
 #include "Curve.h"

 typedef std::vector<ks::vec3> fiber_t;
 typedef std::vector<fiber_t> yarn_t;

 namespace Fiber {

 	/* Define a fiber */
 	struct Fiber {
 		std::vector<ks::vec3> vertices;	 // world-space positions for this fiber
 #ifndef IMPROVED_FLYAWAYS
 		std::vector<std::vector<ks::vec3> > fly_vertices_list; // world-space positions for fly away part
 #endif
 		ks::vec3 init_vertex;				 // initial vertex of this fiber on z=0 normal plane, separate from the previous vertices
 		// because we might need fiber migration 
 		float init_radius;				 // initial radius sampled with Section 4.1 Cross-Sectional Fiber Distribution
 		float init_theta;				 // initial theta for this fiber in circle
 		float init_migration_theta;		 // initial migration theta for this fiber
 		void clear() {
 			vertices.clear();
 #ifndef IMPROVED_FLYAWAYS
 			for (int i = 0; i < fly_vertices_list.size(); i++)
 				fly_vertices_list[i].clear();
 			fly_vertices_list.clear();
 #endif
 		}
 	};

 	/* Define a ply */
 	struct Ply {
 		/* Fibers in this ply */
 		std::vector<Fiber> fibers;

 		std::vector<std::pair<float, float>> cross_section_samples;

 		int num_hair_fibers;

 		/* Base center and theta on z=0 normal plane, and default radius = yarn_radius/2 */
 		ks::vec3 base_center;
 		float base_theta;

 		/* Parameters and functions for Section 4.1 Cross-Sectional Fiber Distribution */
 		float epsilon;
 		float R_max;
 		float beta;
 		float fiberDistribution(float R) {
 			float eTerm = (ks::e - std::powf(ks::e, R / R_max)) / (ks::e - 1);
 			float pR = (1 - 2 * epsilon) * std::powf(eTerm, beta) + epsilon;
 			return pR;
 		}
 		float sampleR() {
 			while (true) {
 				/* [Attention]: Uniformly sample in circle, r = sqrt(rnd1), theta = 2 * pi * rnd2 */
 				/*              This way dS = rdr * dtheta = rnd1^0.5 * 0.5 * rnd1^(-0.5) drnd1 * 2 * pi drnd2 */
 				/*						    = pi * drnd1 * drnd2, which is uniform distributed. */
 				float radius = std::sqrt(ks::rand01()) * R_max, pdf = ks::rand01();
 				if (pdf < fiberDistribution(radius))
 					return radius;
 			}
 		}

 		/* Parameters and functions for Section 4.2 Twist/Fiber Migration */
 		float alpha;
 		float s_i;
 		float rho_min, rho_max;
 		float helixRadius(float init_r, float init_migration_theta, float theta, bool use_migration, float loop_radius) {
 			float r = init_r;

            if (use_migration) {
                if (loop_radius < 0.f) {
                    // triangle wave
                    // float cv = 0;
                    // float migration_theta = init_migration_theta + s_i * theta;
                    // float mt_local = migration_theta - 2 * Zhao::pi * floor(migration_theta / (2 * Zhao::pi));
                    // if (mt_local < Zhao::pi) cv = 1 - 2 * mt_local / Zhao::pi;
                    // else cv = 2 * mt_local / Zhao::pi - 3;
                    // r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f * (cv + 1);
                    r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f *
                        (std::cosf(s_i * (theta)+init_migration_theta) + 1);
                }
                else
                    r = rho_min * init_r +
                    (loop_radius - rho_min * init_r) * 0.5f * (std::cosf(s_i * (theta)+init_migration_theta) + 1);
            }
 			    //r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f * (std::cosf(s_i * (theta + init_migration_theta) ) + 1);
 			//r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f * (std::cosf(s_i * (theta + init_theta) ) + 1);
 		   // r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f * (std::cosf(s_i * (theta)) + 1);
 			return r;
 		}
 		void helixXYZ(float init_r, float init_theta, float theta, bool use_migration, float init_migration_theta, float loop_radius, float& x, float& y/*, float &z*/) {
 			float r = helixRadius(init_r, init_migration_theta, theta, use_migration, loop_radius);
 			x = r * std::cosf(theta + init_theta),
 				y = r * std::sinf(theta + init_theta);
 			/*z = alpha / (2 * pi) * theta;*/
 		}

 		/* Parameters and functions for Section 4.3 Two-Ply Yarns/Strand Compression */
 		float ellipse_long, ellipse_short; // e_i and d_i in original paper

 		bool clock_wise;		 // This parameter determines twisting direction of the fibers in ply 

 #ifdef IMPROVED_FLYAWAYS
 		/* Parameters for Fly-aways */
 		float flyaway_hair_density;
 		float flyaway_hair_ze_mu, flyaway_hair_ze_sigma;
 		float flyaway_hair_r0_mu, flyaway_hair_r0_sigma;
 		float flyaway_hair_re_mu, flyaway_hair_re_sigma;
 		float flyaway_hair_pe_mu, flyaway_hair_pe_sigma;

 		float flyaway_loop_density;
 		float flyaway_loop_r1_mu, flyaway_loop_r1_sigma;
 #else
 		/* Parameters and functions for Section 4.4 Hairiness */
 		float mu, sigma;		 // Flyaway fiber length ~ NormalDistribution(mu, sigma)
 		int flyaway_num;		 // Flyaway fiber number in this ply
 		float fly_step_size;	 // Flyaway step size 
 #endif
 	};

 	class Yarn {
 	public:
 		/* Plys in this yarn */
 		std::vector<Ply> plys;

 		bool use_migration;		// Fiber migration
 		bool clock_wise;		// This parameter determines twisting direction of the plys in yarn 
 		float z_step_size;		// Yarn grows along z-axis, this parameter determines the growth step size
 		float z_step_num;		// Yarn grows with given step number

 		/* Coaxial helix model: x = R(theta)cos(theta), y = R(theta)sin(theta), z = alpha * theta / (2*PI) */
 		float yarn_alpha;		// Plys twisting in yarn
 		float yarn_radius;		// Radius for yarn circle in which plys twist

 		float fiber_radius;

 		bool is_uniform;

 #ifdef IMPROVED_FLYAWAYS
 		bool use_flyaways;
 #endif

 		std::string config_file;  // Config parameters filename
 		std::string output_file;  // Data output filename 

 		Curve z_curve;			  // Fit yarn along the given curve 

 		ks::AABB3 aabb_procedural;	  // Bounding box for all procedural vertices (no flyaway part)
 		ks::AABB3 aabb_micro_ct;		  // Bounding box for Micro-CT volume data

 		omp_lock_t lock;		  // Perlin noise generator needs it

 	public:
 		void setStepSize(const float ss) {
 			this->z_step_size = ss;
 		}
 		void setStepNum(const int sn) {
 			this->z_step_num = (float)sn;
 		}
 		float getStepSize() const {
 			return this->z_step_size;
 		}
 		int getStepNum() const {
 			return (int)this->z_step_num;
 		}
 		int getPlyNum() const {
 			return (int)this->plys.size();
 		}
 		float scaleFacFromRadius(
 			const float newYR, const int newN, const float newSS, const float meshScaler
 		)
 		{
 			std::cout << "scaleFacFromRadius()" << std::endl;
 			float oldYR = this->yarn_radius;
 			float oldSS = this->z_step_size;
 			int oldN = (int)this->z_step_num;
 			float oldLen = oldSS * oldN;

 			float newLen = newSS * newN;

 			float scaleLen = newLen / oldLen;

 			std::cout << "oldLen = " << oldLen << " newLen = " << newLen << " scaleLen=" << scaleLen << std::endl;


 			this->aabb_micro_ct.pMin.z() *= scaleLen;
 			this->aabb_micro_ct.pMax.z() *= scaleLen;

 			const float THE_NEW_SS = 13;

 			this->setStepSize(THE_NEW_SS);
 			this->setStepNum((int)(newLen / THE_NEW_SS));

 			std::cout << "SS=" << this->getStepSize() << " SN=" << this->getStepNum() << std::endl;

 			this->yarn_alpha *= 4; //   twisting
 			this->yarn_alpha *= meshScaler;
 			for (int i = 0; i < (int)plys.size(); i++) {
 				this->plys[i].alpha *= 4; //  twisting
 				this->plys[i].alpha *= meshScaler;
 			}

 			const float scalingFactor = newYR / oldYR;
 			return scalingFactor;

 		}

 	public:
 		Yarn();
 		Yarn(const std::string& path);

 		~Yarn();

 		/* Parse config file from disk */
 		void parse(const char* filename = NULL);

 		/* Save config file to disk*/
 		void save(const char* filename = NULL);

 		/* Simulate yarn */
 		void simulate();
 		/* Simulate ply  */
 		void simulate_ply();
        void simulate_fly_away();

 		/* Load unrolled plys*/
 		ks::vec2 roll_plys(yarn_t& yarn);
 		void roll_plys(const int K, const std::string& ply_fn, const std::string& fiber_fn);

 		/* Write simulated data (separate plys) to disk */
 		void write_plys(const char* filename = NULL);

 		/* Write simulated data (single yarns) to disk */
 		void write_yarn();
 		void write_yarn_binary();
 		void write_binary_fiber_header();
 		void write_total_vertex_num();
 		std::vector<float> fiber_data;

 		/* Simulate woven yarns and write file */
 		int total_fiber_num, total_vertex_num;
 		void reserve_space_for_total_fiber_num();
 		void write_total_fiber_num();
 		//void simulate_and_write_woven(float L1, float L2_h, float L2_w);

 		// kui
 		void createCrossSectionSamples();
 		void createCrossSectionSamplesUniform();
 		void createPoissonSamples();
 		bool conflict(int i, float radius, float theta);

 		std::vector<std::pair<float, float>> crossSectionSamples;
 	};
 }