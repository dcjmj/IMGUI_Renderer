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


#include "fiber.h"



#if 0
#define PERTURB_FIBERS
#ifdef PERTURB_FIBERS
#   define PERTURB_FIBER_PROB          0.9f
#   define PERTURB_FIBER_RATIO         0.25f
#   define PERTURB_FIBER_SMOOTHING     3
#endif
#endif

using namespace ks;

namespace Fiber {
	Yarn::Yarn() {
		this->config_file = "config.txt";
#ifndef INDIVIDUAL_PLY
		this->output_file = "fiber.txt";
#else 
		this->output_file = "individual_ply.txt";
#endif
	}

	Yarn::Yarn(const std::string& path) {
		this->config_file = path + "config.txt";
#ifndef INDIVIDUAL_PLY
		this->output_file = "fiber.txt";
#else 
		this->output_file = "individual_ply.txt";
#endif
	}

	Yarn::~Yarn() {

	}

	/* Treat each ply with the same settings */
	void Yarn::parse(const char* filename) {
		std::ifstream fin;
		if (filename != NULL)
			fin.open(filename);
		else
			fin.open(this->config_file.c_str());

		std::string line;
		while (std::getline(fin, line)) {
#ifdef VERBOSE
			std::cout << line << std::endl;
#endif
			std::vector<std::string> splits = split(line, ' ');
			if (splits.size() < 2)    continue;
			std::string p_name = splits[0];

			if (p_name == "ply_num:") {
				this->plys.resize(atoi(splits[1].c_str()));
			}
			else if (p_name == "fiber_num:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					int fiber_num = atoi(splits[1].c_str());
					this->plys[i].fibers.resize(fiber_num);
				}
#ifdef IMPROVED_FLYAWAYS
			}
			else if (p_name == "use_flyaways:") {
				this->use_flyaways = (atoi(splits[1].c_str()) != 0);
			}
			else if (p_name == "flyaway_hair_density:") {
				assert(this->plys.size());
				assert(splits.size() == 2);
				for (int i = 0; i < (int)this->plys.size(); i++)
					this->plys[i].flyaway_hair_density = (float)atof(splits[1].c_str());
			}
			else if (p_name == "flyaway_hair_ze:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < (int)this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_ze_mu = (float)atof(splits[1].c_str());
					this->plys[i].flyaway_hair_ze_sigma = (float)atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_r0:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < (int)this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_r0_mu = (float)atof(splits[1].c_str());
					this->plys[i].flyaway_hair_r0_sigma = (float)atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_re:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < (int)this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_re_mu = (float)atof(splits[1].c_str());
					this->plys[i].flyaway_hair_re_sigma = (float)atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_pe:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < (int)this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_pe_mu = (float)atof(splits[1].c_str());
					this->plys[i].flyaway_hair_pe_sigma = (float)atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_loop_density:") {
				assert(this->plys.size());
				assert(splits.size() == 2);
				for (int i = 0; i < (int)this->plys.size(); i++)
					this->plys[i].flyaway_loop_density = (float)atof(splits[1].c_str());
			}
			else if (p_name == "flyaway_loop_r1:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < (int)this->plys.size(); i++)
				{
					this->plys[i].flyaway_loop_r1_mu = (float)atof(splits[1].c_str());
					this->plys[i].flyaway_loop_r1_sigma = (float)atof(splits[2].c_str());
				}
#else
			}
			else if (p_name == "flyaway_num:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].flyaway_num = atoi(splits[1].c_str());
			}
			else if (p_name == "fly_step_size:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].fly_step_size = atof(splits[1].c_str());
#endif
			}
			else if (p_name == "z_step_size:") {
				this->z_step_size = (float)atof(splits[1].c_str());
			}
			else if (p_name == "z_step_num:") {
				this->z_step_num = (float)atof(splits[1].c_str());
			}
			else if (p_name == "yarn_clock_wise:") {
				this->clock_wise = (atoi(splits[1].c_str()) != 0);
			}
			else if (p_name == "fiber_clock_wise:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++)
					this->plys[i].clock_wise = (atoi(splits[1].c_str()) != 0);
			}
			else if (p_name == "use_migration:") {
				this->use_migration = (atoi(splits[1].c_str()) != 0);
			}
			else if (p_name == "yarn_alpha:") {
				this->yarn_alpha = (float)atof(splits[1].c_str());
			}
			else if (p_name == "yarn_radius:") {
				this->yarn_radius = (float)atof(splits[1].c_str());
			}
			else if (p_name == "epsilon:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].epsilon = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "R_max:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].R_max = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "beta:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].beta = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "alpha:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].alpha = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "s_i:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].s_i = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "rho_min:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].rho_min = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "rho_max:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].rho_max = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "ellipse_long:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].ellipse_long = (float)atof(splits[1].c_str());
				}
			}
			else if (p_name == "ellipse_short:") {
				assert(this->plys.size());
				for (int i = 0; i < (int)this->plys.size(); i++) {
					this->plys[i].ellipse_short = (float)atof(splits[1].c_str());
				}
#ifndef IMPROVED_FLYAWAYS
			}
			else if (p_name == "mu:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].mu = atof(splits[1].c_str());
				}
			}
			else if (p_name == "sigma:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].sigma = atof(splits[1].c_str());
				}
#endif
			}
			else if (p_name == "z_curve_file:") {
				//this->z_curve.setFile(splits[1]); // Only one curve
			}
			else if (p_name == "z_curves_file:") {
				//this->z_curve.setFile(splits[1]);
				//this->z_curve.setMultiCurveFlag(true); // Multiple curves (e.g., woven)
			}
			else if (p_name == "aabb_min:") {
				std::string min_str = splits[1];
				std::vector<std::string> min_values = split(min_str.substr(1, min_str.size() - 2), ',');
				assert(min_values.size() == 3);
				for (int i = 0; i < 3; i++) {
					this->aabb_micro_ct.pMin[i] = (float)atof(min_values[i].c_str());
				}
			}
			else if (p_name == "aabb_max:") {
				std::string max_str = splits[1];
				std::vector<std::string> max_values = split(max_str.substr(1, max_str.size() - 2), ',');
				assert(max_values.size() == 3);
				for (int i = 0; i < 3; i++) {
					this->aabb_micro_ct.pMax[i] = (float)atof(max_values[i].c_str());
				}
			}
		}
		fin.close();
	}

	void Yarn::save(const char* filename) {
		std::ofstream fout;
		if (filename != NULL)
			fout.open(filename);
		else
			fout.open(this->config_file.c_str());

		fout << "ply_num: " << this->plys.size() << std::endl;
		fout << "fiber_num: " << this->plys[0].fibers.size() << std::endl;
#ifndef IMPROVED_FLYAWAYS
		fout << "flyaway_num: " << this->flyaway_num << std::endl;
#endif
		fout << std::endl;

		fout << "aabb_min: [" << this->aabb_micro_ct.pMin.x() << "," << this->aabb_micro_ct.pMin.y() << "," << this->aabb_micro_ct.pMin.z() << "]" << std::endl;
		fout << "aabb_max: [" << this->aabb_micro_ct.pMax.x() << "," << this->aabb_micro_ct.pMax.y() << "," << this->aabb_micro_ct.pMax.z() << "]" << std::endl;
		fout << std::endl;

		fout << "z_step_size: " << this->z_step_size << std::endl;
		fout << "z_step_num: " << this->z_step_num << std::endl;
#ifndef IMPROVED_FLYAWAYS
		fout << "fly_step_size: " << this->fly_step_size << std::endl;
#endif
		fout << std::endl;

		fout << "yarn_clock_wise: " << (this->clock_wise ? 1 : 0) << std::endl;
		fout << "fiber_clock_wise: " << (this->plys[0].clock_wise ? 1 : 0) << std::endl;
		fout << "yarn_alpha: " << this->yarn_alpha << std::endl;
		fout << "alpha: " << this->plys[0].alpha << std::endl;
		fout << std::endl;

		fout << "yarn_radius: " << this->yarn_radius << std::endl;
		fout << "ellipse_long: " << this->plys[0].ellipse_long << std::endl;
		fout << "ellipse_short: " << this->plys[0].ellipse_short << std::endl;
		fout << std::endl;

		fout << "epsilon: " << this->plys[0].epsilon << std::endl;
		fout << "beta: " << this->plys[0].beta << std::endl;
		fout << "R_max: " << this->plys[0].R_max << std::endl;
		fout << std::endl;

		fout << "use_migration: " << this->use_migration << std::endl;
		fout << "s_i: " << this->plys[0].s_i << std::endl;
		fout << "rho_min: " << this->plys[0].rho_min << std::endl;
		fout << "rho_max: " << this->plys[0].rho_max << std::endl;
		fout << std::endl;

#ifdef IMPROVED_FLYAWAYS
		fout << "use_flyaways: " << this->use_flyaways << std::endl;
		fout << "flyaway_hair_density: " << this->plys[0].flyaway_hair_density << std::endl;
		fout << "flyaway_hair_ze: " << this->plys[0].flyaway_hair_ze_mu << ' ' << this->plys[0].flyaway_hair_ze_sigma << std::endl;
		fout << "flyaway_hair_r0: " << this->plys[0].flyaway_hair_r0_mu << ' ' << this->plys[0].flyaway_hair_r0_sigma << std::endl;
		fout << "flyaway_hair_re: " << this->plys[0].flyaway_hair_re_mu << ' ' << this->plys[0].flyaway_hair_re_sigma << std::endl;
		fout << "flyaway_hair_pe: " << this->plys[0].flyaway_hair_pe_mu << ' ' << this->plys[0].flyaway_hair_pe_sigma << std::endl;
		fout << "flyaway_loop_density: " << this->plys[0].flyaway_loop_density << std::endl;
		fout << "flyaway_loop_r1: " << this->plys[0].flyaway_loop_r1_mu << ' ' << this->plys[0].flyaway_loop_r1_sigma << std::endl;
#else
		fout << "mu: " << this->mu << std::endl;
		fout << "sigma: " << this->sigma << std::endl;
#endif

		fout.close();
	}

	void Yarn::simulate() {
		omp_init_lock(&this->lock);
		// Step1: Obtain center of yarn starting point 
		printf("Obtain center of yarn starting point...\n");
		const vec3 base_center = vec3(0, 0, 0);
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();
		// Step2: Sample initial locations of ply-centers in normal plane around starting point
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		// Step3: Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1
		printf("Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1...\n");
		const int num_of_cores = omp_get_num_procs();

		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];
				float radius = this->plys[i].sampleR();
				float theta = 2 * pi * rand01();
				float migration_theta = 2 * pi * rand01();
				fiber.init_radius = radius;
				fiber.init_theta = theta;
				fiber.init_migration_theta = migration_theta;
				fiber.init_vertex = this->plys[i].base_center +
					vec3(radius * std::cosf(theta), radius * std::sinf(theta), 0);
			}

		}

		//// Added: Fit curve along given curve, otherwise along z-axis
		//if (this->z_curve.hasCurve()) {
		//	this->z_curve.load();
		//	this->z_curve.setStepNumbers(this->z_step_num);
		//}

		// Step4: Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers...\n");
#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];
				fiber.clear();

#ifdef PERTURB_FIBERS
				std::vector<float> perturbRatios;
				{
					const int nsteps = static_cast<int>(std::ceil(this->z_step_num));
					std::vector<int> eventLoc;
					for (int step_id = 0; step_id < nsteps; ++step_id)
						if (rand01() < PERTURB_FIBER_PROB)
							eventLoc.push_back(step_id);
					perturbRatios.resize(nsteps, 1.0f);
					if (!eventLoc.empty())
					{
						std::vector<int>::iterator it = eventLoc.begin();
						perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO * (rand01() - 0.5f);
						for (int j = 0; j < *it; ++j) perturbRatios[j] = perturbRatios[*it];
						while ((++it) != eventLoc.end())
						{
							perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO * (rand01() - 0.5f);
							float extent = static_cast<float>(*it - *(it - 1));
							for (int j = *(it - 1) + 1; j < *it; ++j)
								perturbRatios[j] = (perturbRatios[*(it - 1)] * (*it - j) + perturbRatios[*it] * (j - *(it - 1))) / extent;
						}
						for (int j = eventLoc.back() + 1; j < nsteps; ++j)
							perturbRatios[j] = perturbRatios[eventLoc.back()];
					}

					for (int j = 0; j < PERTURB_FIBER_SMOOTHING; ++j)
					{
						std::vector<float> perturbRatios0 = perturbRatios;
						for (int k = 1; k + 1 < nsteps; ++k)
							perturbRatios[k] = 0.25f * perturbRatios0[k - 1] + 0.5f * perturbRatios0[k] + 0.25f * perturbRatios0[k + 1];
					}
				}
#endif

				for (int step_id = 0; step_id < this->z_step_num; step_id++) {
					const float z = this->z_step_size * (step_id - this->z_step_num / 2.f);
					const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float yarn_theta = this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
					float local_x, local_y, world_x, world_y;

					// Step5: Vary the distance of cross-sectional fiber positions to their ply center according to fiber migration Sec 4.2
					this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, -1, local_x, local_y);
#ifndef INDIVIDUAL_PLY
					// Step 6: Transform cross-sectional fiber positions according to strand compression Sec 4.3
					vec3 short_axis = (this->plys[i].base_center).normalized(), long_axis = vec3(-short_axis.y(), short_axis.x(), 0);
					vec3 local_p = vec3(local_x, local_y, 0.f);
					float _local_x = local_p.dot(short_axis), _local_y = local_p.dot(long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x();
					local_y = local_p.y();

#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x();
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y();
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);
#else 
					const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);
					local_x *= balance_radius;
					local_y *= balance_radius;
#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif
					world_x = local_x;
					world_y = local_y;
#endif
					vec3 verIn = vec3(world_x, world_y, z), verOut;

					//// Added: Fit curve along given curve, otherwise along z-axis
					//if (this->z_curve.hasCurve()) {
					//	this->z_curve.fit(step_id, verIn, verOut);
					//} else {
					verOut = verIn;
					//}

					this->aabb_procedural.expand(verOut);
					if (this->aabb_micro_ct.contain(verOut))
						fiber.vertices.push_back(verOut);
				}
			}
		}

		std::cout << "Bounding box:\n"
			<< "aabb_min = [" << this->aabb_procedural.pMin.x() << "," << this->aabb_procedural.pMin.y() << "," << this->aabb_procedural.pMin.z() << "]" << std::endl
			<< "aabb_max = [" << this->aabb_procedural.pMax.x() << "," << this->aabb_procedural.pMax.y() << "," << this->aabb_procedural.pMax.z() << "]" << std::endl;

#ifndef IMPROVED_FLYAWAYS
		// Step 7: Transform cross-sectional fiber positions according to strand compression Sec 4.3
		printf("Choose flyaway fibers as in Sec 4.4...\n");
#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {
			const int fiber_flyaway_num = this->plys[i].flyaway_num;
			this->plys[i].fly_fiber_num = 0;
			std::default_random_engine generator;
			std::normal_distribution<float> distribution(this->plys[i].mu, this->plys[i].sigma);
			PerlinNoise noise1, noise2, noise3;
			for (int f = 0; f < fiber_flyaway_num; ) {
				/* Randomly uniformly pick one vertex from fibers */
				const int fiber_chosen = std::floor(this->plys[i].fibers.size() * rand01());
				Fiber& fiber = this->plys[i].fibers[fiber_chosen];

				const int vertex_chosen = std::floor(fiber.vertices.size() * rand01());
				if (vertex_chosen == 0 || vertex_chosen == fiber.vertices.size() - 1)
					continue;
				const vec3 vertex = fiber.vertices[vertex_chosen];
				const int vertex_chosen_next = vertex_chosen + 1;
				const vec3 vertex_next = fiber.vertices[vertex_chosen_next];
				/* Find world-space tangent direction */
				bool coin_flip_result = coin_flip();
				vec3 tangent = coin_flip_result ?
					(vertex_next - vertex).normalized() : (vertex - vertex_next).normalized();

				/* Sample a distance according to the normal distribution */
				const float distance = distribution(generator);
				const int num_steps = (int)std::max(0.f, std::ceil(distance / this->plys[i].fly_step_size));

				if (num_steps > 0) {
					std::vector<vec3> fly_vertices;
					const float flyaway_step_size = distance / num_steps;

					vec3 flow_vertex = this->aabb_procedural.shift(vertex, flyaway_step_size);

					int step_id = 0;
					while (step_id++ < num_steps) {
						fly_vertices.push_back(flow_vertex);
						omp_set_lock(&lock);
						vec3 perlin = perlin_vector_field(
							noise1, noise2, noise3,
							flow_vertex, this->yarn_radius
						);
						omp_unset_lock(&lock);
						flow_vertex += (perlin + tangent) * flyaway_step_size;
						if (/*this->aabb_procedural.in(flow_vertex) || */this->aabb_micro_ct.out(flow_vertex))
							break;
					}

					omp_set_lock(&lock);
					if (fly_vertices.size() >= 2) {
						fiber.fly_vertices_list.push_back(fly_vertices);
						this->plys[i].fly_fiber_num++;
						f++;
					}
					omp_unset_lock(&lock);
				}
			}
		}
		omp_destroy_lock(&this->lock);
#endif
		printf("Simulation is done!\n");
	}

	void Yarn::simulate_ply() {

		std::cout << "simulate_ply()" << std::endl;

#define INDIVIDUAL_PLY 
		omp_init_lock(&this->lock);
		// Step1: Obtain center of yarn starting point 
#ifdef VERBOSE
		printf("Obtain center of yarn starting point...\n");
#endif
		const vec3 base_center = vec3(0, 0, 0);
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();
		// Step2: Sample initial locations of ply-centers in normal plane around starting point
#ifdef VERBOSE
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
#endif
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		// Step3: Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1
#ifdef VERBOSE
		printf("Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1...\n");
#endif

		const int num_of_cores = omp_get_num_procs();
		srand(0);
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
//#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];

				// TODO:
				float radius, theta;
				if (is_uniform)
				{
					radius = 1.0f;
					theta = 2 * pi * f / fiber_num;
				}
				else
				{
					radius = this->plys[i].sampleR();
					theta = 2 * pi * rand01();
				}

				float migration_theta = 2 * pi * rand01();
				fiber.init_radius = radius;
				fiber.init_theta = theta;
				fiber.init_migration_theta = migration_theta;
				fiber.init_vertex = this->plys[i].base_center +
					vec3(radius * std::cosf(theta), radius * std::sinf(theta), 0);

				std::pair<float, float> p(radius, theta);
				this->plys[i].cross_section_samples.push_back(p);
				//std::cout << i << " " << f << " " << theta << " " << radius << std::endl;
			}

		}


		// Step4: Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers
#ifdef VERBOSE	
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers...\n");
#endif

#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();

			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];
				fiber.clear();

#ifdef PERTURB_FIBERS
				std::vector<float> perturbRatios;
				{
					const int nsteps = static_cast<int>(std::ceil(this->z_step_num));
					std::vector<int> eventLoc;
					for (int step_id = 0; step_id < nsteps; ++step_id)
						if (rand01() < PERTURB_FIBER_PROB)
							eventLoc.push_back(step_id);
					perturbRatios.resize(nsteps, 1.0f);
					if (!eventLoc.empty())
					{
						std::vector<int>::iterator it = eventLoc.begin();
						perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO * (rand01() - 0.5f);
						for (int j = 0; j < *it; ++j) perturbRatios[j] = perturbRatios[*it];
						while ((++it) != eventLoc.end())
						{
							perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO * (rand01() - 0.5f);
							float extent = static_cast<float>(*it - *(it - 1));
							for (int j = *(it - 1) + 1; j < *it; ++j)
							{
#if 0
								perturbRatios[j] = (perturbRatios[*(it - 1)] * (*it - j) + perturbRatios[*it] * (j - *(it - 1))) / extent;
#else
								float v = static_cast<float>(*it - j) / extent;
								v = std::sin(0.5f * pi * v);
								perturbRatios[j] = perturbRatios[*(it - 1)] * v + perturbRatios[*it] * (1.0f - v);
#endif
							}
						}
						for (int j = eventLoc.back() + 1; j < nsteps; ++j)
							perturbRatios[j] = perturbRatios[eventLoc.back()];
					}

					for (int j = 0; j < PERTURB_FIBER_SMOOTHING; ++j)
					{
						std::vector<float> perturbRatios0 = perturbRatios;
						for (int k = 1; k + 1 < nsteps; ++k)
							perturbRatios[k] = 0.25f * perturbRatios0[k - 1] + 0.5f * perturbRatios0[k] + 0.25f * perturbRatios0[k + 1];
					}
				}
#endif

				for (int step_id = 0; step_id < this->z_step_num - 1; step_id++) {
					float z = this->z_step_size * (step_id - this->z_step_num / 2.f);
					float z_real = z;
					if (step_id > this->z_step_num - 2) {
						z_real = this->z_step_size * (this->z_step_num - 2) / 2.f;
						z = this->z_step_size * (this->z_step_num - 2) / 2.f;
					}
					const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float yarn_theta = this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
					float local_x, local_y, world_x, world_y;


					// Step5: Vary the distance of cross-sectional fiber positions to their ply center according to fiber migration Sec 4.2
					this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, -1, local_x, local_y);
#ifndef INDIVIDUAL_PLY
					// Step 6: Transform cross-sectional fiber positions according to strand compression Sec 4.3
					vec3 short_axis = (this->plys[i].base_center).normalized(), long_axis = vec3(-short_axis.y, short_axis.x, 0);
					vec3 local_p = vec3(local_x, local_y, 0.f);
					float _local_x = local_p.dot(short_axis), _local_y = local_p.dot(long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x;
					local_y = local_p.y;

#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x;
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y;
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);
#else 
					const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);
					local_x *= balance_radius;
					local_y *= balance_radius;

#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					world_x = local_x;
					world_y = local_y;
#endif
					vec3 verIn = vec3(world_x, world_y, z_real);
					fiber.vertices.push_back(verIn);
					/*if ((step_id == 0) || (step_id > this->z_step_num - 2)) {
						std::cout << step_id << " " << fiber_theta << " " << yarn_theta << std::endl;
						std::cout << verIn.x() << " " << verIn.y() << " " << verIn.z() << std::endl;
					}*/
				}
			}
		}

		for (int i = 0; i < ply_num; ++i)
			for (auto it = this->plys[i].fibers.begin(); it != this->plys[i].fibers.end(); ++it)
				for (auto it2 = it->vertices.begin(); it2 != it->vertices.end(); ++it2) {
					this->aabb_procedural.expand(*it2);
				}


		omp_destroy_lock(&this->lock);

		std::cout << "Checking..." << std::endl;
		int bad_count = 0;
		for (int i = 0; i < (int)this->plys.size(); i++)
		{
			const int fiberNum = this->plys[i].fibers.size();
			for (int f = 0; f < fiberNum; f++)
			{
				const int vertexNum = this->plys[i].fibers[f].vertices.size();
				for (int v = 1; v < vertexNum - 1; v++)
				{
					vec3 prev = this->plys[i].fibers[f].vertices[v - 1];
					vec3 curr = this->plys[i].fibers[f].vertices[v];
					vec3 next = this->plys[i].fibers[f].vertices[v + 1];
					vec3 dir1 = (next - curr).normalized();
					vec3 dir2 = (curr - prev).normalized();
					if (dir1.dot(dir2) < 0.5) {
						bad_count++;
					}
				}
			}
		}
		std::cout << "Bad count = " << bad_count << std::endl;

	}

	vec2 Yarn::roll_plys(yarn_t& r_yarn) {
		std::cout << "roll_plys() " << std::endl;

		const int num_of_cores = omp_get_num_procs();
		const vec3 base_center = vec3(0, 0, 0);
		const float base_radius = this->yarn_radius;

		const int ply_num = this->plys.size();

		int total_fiber_num = 0;
		std::vector<int> accumulate_ply_fibers_num(ply_num);


		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);

			accumulate_ply_fibers_num[i] = total_fiber_num;
			total_fiber_num += this->plys[i].fibers.size();
		}

		r_yarn.resize(total_fiber_num);


		const float yarn_alpha = this->yarn_alpha;
		const int   yarn_clock = this->clock_wise;

		vec2 z_minmax = vec2(FLT_MAX, FLT_MIN);
		omp_lock_t lock;
		omp_init_lock(&lock);

		for (int i = 0; i < ply_num; i++) {
			const float e_l = this->plys[i].ellipse_long, e_s = this->plys[i].ellipse_short, b_r = sqrtf(e_l * e_s);
			const int fiberNum = this->plys[i].fibers.size();
			for (int f = 0; f < fiberNum; f++) {
				int fid = accumulate_ply_fibers_num[i] + f;

				const int vertexNum = this->plys[i].fibers[f].vertices.size();
				r_yarn[fid].resize(vertexNum);
#pragma omp parallel for num_threads(num_of_cores) 
				for (int v = 0; v < vertexNum; v++) {
					const float z = this->plys[i].fibers[f].vertices[v].z();     //yarn[i][f][v].z;
					const float yarn_theta = yarn_clock ? -z * 2 * pi / yarn_alpha : z * 2 * pi / yarn_alpha;
					float local_x, local_y, world_x, world_y;

					local_x = this->plys[i].fibers[f].vertices[v].x() / b_r;
					local_y = this->plys[i].fibers[f].vertices[v].y() / b_r;


					vec3 short_axis = (this->plys[i].base_center).normalized(), long_axis = vec3(-short_axis.y(), short_axis.x(), 0);
					vec3 local_p = vec3(local_x, local_y, 0.f);
					float _local_x = local_p.dot(short_axis), _local_y = local_p.dot(long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x();
					local_y = local_p.y();
					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x();
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y();
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);

					this->plys[i].fibers[f].vertices[v].x() = world_x;
					this->plys[i].fibers[f].vertices[v].y() = world_y;
					r_yarn[fid][v] = vec3(world_x, world_y, z);

					omp_set_lock(&lock);
					z_minmax.x() = std::min(z_minmax.x(), z);
					z_minmax.y() = std::max(z_minmax.y(), z);
					omp_unset_lock(&lock);
				}
			}

		}

		std::ofstream fout("fiber_before_warping.txt");
		fiber_t& the_fiber = r_yarn[0];
		int the_fiber_vertexNum = the_fiber.size();
		for (int v = 0; v < the_fiber_vertexNum; v++)
		{
			fout << the_fiber[v].x() << " " << the_fiber[v].y() << " " << the_fiber[v].z() << std::endl;
		}

		fout.close();


		std::cout << "Checking..." << std::endl;
		int bad_count = 0;
		const int fNum = r_yarn.size();
		for (int f = 0; f < fNum; f++)
		{
			fiber_t& fiber = r_yarn[f];
			const int vNum = fiber.size();
			for (int v = 1; v < vNum - 1; v++)
			{
				vec3 prev = fiber[v - 1];
				vec3 curr = fiber[v];
				vec3 next = fiber[v + 1];
				vec3 dir1 = (next - curr).normalized();
				vec3 dir2 = (curr - prev).normalized();
				if (dir1.dot(dir2) < 0.5) {
					bad_count++;
				}
			}
		}
		std::cout << "Bad count = " << bad_count << std::endl;
		omp_destroy_lock(&lock);
		return z_minmax;

	}

	void Yarn::roll_plys(const int K, const std::string& ply_fn, const std::string& fiber_fn) {
		const int num_of_cores = omp_get_num_procs();
#ifdef VERBOSE
		std::cout << "Rolling plys into yarn... " << std::endl;
#endif
		typedef std::vector<vec3> fiber_t;
		typedef std::vector<fiber_t> ply_t;

		std::vector<ply_t> yarn(K);
		int total_fiber_num = 0; std::vector<int> ply_fiber_num(K);

		for (int i = 0; i < K; i++) {
			std::string filename = ply_fn.substr(0, ply_fn.find('.')) + std::to_string((long long)i) + ".txt";
#ifdef VERBOSE
			std::cout << "Loading unrolled ply = " << filename << std::endl;
#endif
			std::ifstream fin(filename.c_str());
			int fiberNum; fin >> fiberNum;
			ply_fiber_num[i] = fiberNum;
			total_fiber_num += fiberNum;
			yarn[i].resize(fiberNum);
			for (int f = 0; f < fiberNum; f++) {
				int vertexNum; fin >> vertexNum;
				yarn[i][f].resize(vertexNum);
				for (int v = 0; v < vertexNum; v++) {
					fin >> yarn[i][f][v].x() >> yarn[i][f][v].y() >> yarn[i][f][v].z();
				}
			}
			fin.close();
		}

#ifdef VERBOSE
		printf("Obtain center of yarn starting point...\n");
#endif
		const vec3 base_center = vec3(0, 0, 0);
		const float base_radius = this->yarn_radius;

#ifdef VERBOSE
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
#endif
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		const float yarn_alpha = this->yarn_alpha;
		const int   yarn_clock = this->clock_wise;

#ifdef VERBOSE
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center...\n");
#endif
#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {

			const float e_l = this->plys[i].ellipse_long, e_s = this->plys[i].ellipse_short, b_r = sqrtf(e_l * e_s);

			const int fiberNum = yarn[i].size();
			for (int f = 0; f < fiberNum; f++) {
				const int vertexNum = yarn[i][f].size();
				for (int v = 0; v < vertexNum; v++) {
					const float z = yarn[i][f][v].z();
					const float yarn_theta = yarn_clock ? -z * 2 * pi / yarn_alpha : z * 2 * pi / yarn_alpha;
					float local_x, local_y, world_x, world_y;

					local_x = yarn[i][f][v].x() / b_r;
					local_y = yarn[i][f][v].y() / b_r;

					//std::cout << this->plys[i].base_center.x << " " << this->plys[i].base_center.y << " " << this->plys[i].base_center.z << std::endl;
					vec3 short_axis = (this->plys[i].base_center).normalized(), long_axis = vec3(-short_axis.y(), short_axis.x(), 0);
					vec3 local_p = vec3(local_x, local_y, 0.f);
					float _local_x = local_p.dot(short_axis), _local_y = local_p.dot(long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x();
					local_y = local_p.y();
					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x();
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y();
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);

					yarn[i][f][v].x() = world_x;
					yarn[i][f][v].y() = world_y;
				}
			}

		}

#ifdef VERBOSE
		std::cout << "Writing final yarn file... " << std::endl;
#endif
		std::ofstream fout(fiber_fn.c_str());

		fout << total_fiber_num << std::endl;

		for (int i = 0; i < K; i++) {
			const int fiberNum = yarn[i].size();
			for (int f = 0; f < fiberNum; f++) {
				int vertexNum = yarn[i][f].size();
				fout << vertexNum << std::endl;
				for (int v = 0; v < vertexNum; v++) {
					fout << yarn[i][f][v].x() << " " << yarn[i][f][v].y() << " " << yarn[i][f][v].z() << std::endl;
				}
			}
		}


		fout.close();
#ifdef VERBOSE
		std::cout << "Done!" << std::endl;
#endif
	}

	void Yarn::write_plys(const char* fn) {
#ifdef VERBOSE
		printf("Writing vertices to file...\n");
#endif
		std::string filename;
		int ply_num = this->plys.size();
		std::string path = WORK_PATH;
		std::vector<std::ofstream> fouts(ply_num);
		for (int i = 0; i < ply_num; i++) {
			if (fn) {
				filename = std::string(fn);
				filename = filename.substr(0, filename.find('.')) + std::to_string(static_cast<long long>(i)) + ".txt";
			}
			else
				filename = (path + this->output_file.substr(0, this->output_file.find('.')) + std::to_string(static_cast<long long>(i)) + ".txt");

			fouts[i].open(filename.c_str());
#ifdef VERBOSE
			std::cout << "Writing File : " << filename << std::endl;
#endif
			int fiber_num = this->plys[i].fibers.size();
#ifndef IMPROVED_FLYAWAYS
			int fly_fiber_num = this->plys[i].fly_fiber_num;
#endif

#ifdef VERBOSE
			std::cout << "fiber_num = " << fiber_num
#ifndef IMPROVED_FLYAWAYS
				<< " flyaway_num = " << fly_fiber_num
#endif
				<< std::endl;
#endif

			fouts[i] << fiber_num
#ifndef IMPROVED_FLYAWAYS
				+ fly_fiber_num
#endif
				<< std::endl;

			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];

				int fiber_vertex_num = fiber.vertices.size();
				fouts[i] << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fouts[i] << fiber.vertices[v].x() << " " << fiber.vertices[v].y() << " " << fiber.vertices[v].z() << std::endl;
				}

#ifndef IMPROVED_FLYAWAYS
				int fly_fiber_num = fiber.fly_vertices_list.size();
				if (fly_fiber_num > 0) {
					for (int fi = 0; fi < fly_fiber_num; fi++) {
						int fly_fiber_vertex = fiber.fly_vertices_list[fi].size();
						fouts[i] << fly_fiber_vertex << std::endl;
						for (int v = 0; v < fly_fiber_vertex; v++)
							fouts[i] << fiber.fly_vertices_list[fi][v].x << " " << fiber.fly_vertices_list[fi][v].y << " " << fiber.fly_vertices_list[fi][v].z << std::endl;
					}
				}
#endif
			}
			fouts[i].close();
		}

#ifdef VERBOSE
		printf("Writing vertices to file done!\n");
#endif
	}

	void Yarn::write_yarn() {
		bool append = false; //this->z_curve.multiCurves();
		printf("Writing vertices to file...\n");
		std::string path = WORK_PATH;
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
#ifdef IMPROVED_FLYAWAYS
			total_fiber_num += this->plys[i].fibers.size();
#else
			total_fiber_num += this->plys[i].fibers.size() + this->plys[i].fly_fiber_num;
#endif
		std::ofstream fout((path + this->output_file).c_str(),
			append ? std::ofstream::app : std::ofstream::trunc);
		if (append)
			this->total_fiber_num += total_fiber_num;
		else
			fout << total_fiber_num << std::endl;

		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				fout << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fout << fiber.vertices[v].x() << " " << fiber.vertices[v].y() << " " << fiber.vertices[v].z() << std::endl;
				}

#ifndef IMPROVED_FLYAWAYS
				int fly_fiber_num = fiber.fly_vertices_list.size();
				if (fly_fiber_num > 0) {
					for (int fi = 0; fi < fly_fiber_num; fi++) {
						int fly_fiber_vertex = fiber.fly_vertices_list[fi].size();
						fout << fly_fiber_vertex << std::endl;
						for (int v = 0; v < fly_fiber_vertex; v++)
							fout << fiber.fly_vertices_list[fi][v].x << " " << fiber.fly_vertices_list[fi][v].y << " " << fiber.fly_vertices_list[fi][v].z << std::endl;
					}
				}
#endif
			}
		}
		fout.close();
		printf("Writing vertices to file done!\n");
	}

	void Yarn::write_yarn_binary() {
		bool append = false;//this->z_curve.multiCurves();
		printf("Writing vertices to file...\n");
		std::string path = WORK_PATH;
		fiber_data.clear();
		int total_vertex_num = 0, total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
#ifndef IMPROVED_FLYAWAYS
			int fly_fiber_num = this->plys[i].fly_fiber_num;
			total_fiber_num += fiber_num + fly_fiber_num;
#endif
			for (int f = 0; f < fiber_num; f++) {
				Fiber& fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				fiber_data.push_back(std::numeric_limits<float>::infinity());
				for (int v = 0; v < fiber_vertex_num; v++) {
					fiber_data.push_back(fiber.vertices[v].x());
					fiber_data.push_back(fiber.vertices[v].y());
					fiber_data.push_back(fiber.vertices[v].z());
				}
				total_vertex_num += fiber_vertex_num;

#ifndef IMPROVED_FLYAWAYS
				int fly_fiber_num = fiber.fly_vertices_list.size();
				if (fly_fiber_num > 0) {
					for (int fi = 0; fi < fly_fiber_num; fi++) {
						int fly_fiber_vertex = fiber.fly_vertices_list[fi].size();
						total_vertex_num += fly_fiber_vertex;
						for (int v = 0; v < fly_fiber_vertex; v++) {
							fiber_data.push_back(fiber.fly_vertices_list[fi][v].x());
							fiber_data.push_back(fiber.fly_vertices_list[fi][v].y());
							fiber_data.push_back(fiber.fly_vertices_list[fi][v].z());
						}
					}
				}
#endif
			}
		}
		this->total_fiber_num += total_fiber_num;
		this->total_vertex_num += total_vertex_num;
		FILE* file = fopen((path + this->output_file).c_str(), "ab");
		fwrite(&fiber_data[0], sizeof(float), fiber_data.size(), file);
		fclose(file);
		printf("Writing vertices to file done!\n");
	}

	void Yarn::reserve_space_for_total_fiber_num() {
		std::string path = WORK_PATH;
		const std::string reserve_string = "        ";
		std::ofstream fout((path + this->output_file).c_str());
		fout << reserve_string << std::endl;
		fout.close();
	}

	void Yarn::write_total_fiber_num() {
		std::string path = WORK_PATH;
		std::fstream file((path + this->output_file).c_str(), std::ios::in | std::ios::out);
		file.seekp(0);
		file << this->total_fiber_num;
		file.close();
	}

	void Yarn::write_total_vertex_num() {
		std::string path = WORK_PATH;
		const int header_len = 11;
		FILE* file = fopen((path + this->output_file).c_str(), "r+");
		fseek(file, header_len, SEEK_SET);
		fwrite(&total_vertex_num, sizeof(int), 1, file);
		fclose(file);
	}

	void Yarn::write_binary_fiber_header() {
		std::string path = WORK_PATH;
		const int header_len = 11;
		const char header[12] = "BINARY_HAIR";
		FILE* file = fopen((path + this->output_file).c_str(), "wb");
		/* Mitsuba binary fiber header */
		fwrite(header, sizeof(char), header_len, file);
		/* Reserve total vertex number int */
		int total_vertex_num = 0;
		fwrite(&total_vertex_num, sizeof(int), 1, file);
		fclose(file);
	}
}

bool Fiber::Yarn::conflict(int i, float radius, float theta)
{
	for (int j = 0; j < (int)this->plys[i].cross_section_samples.size(); j++)
	{
		vec2 p(cos(theta) * radius, sin(theta) * radius);
		vec2 tmpp(cos(this->plys[i].cross_section_samples[j].second) * this->plys[i].cross_section_samples[j].first,
			sin(this->plys[i].cross_section_samples[j].second) * this->plys[i].cross_section_samples[j].first);
		if ((p - tmpp).norm() < fiber_radius)
		{
			return true;
		}
	}
	return false;
}

void Fiber::Yarn::createPoissonSamples()
{
	std::cout << "Create Poisson Samples\n";
	float radius, theta;
	const int ply_num = this->plys.size();
	for (int j = 0; j < ply_num; j++) {
		this->plys[j].cross_section_samples.clear();
		radius = this->plys[0].sampleR();
		theta = 2 * pi * rand01();
		std::pair<float, float> sample(radius, theta);
		this->plys[j].cross_section_samples.push_back(sample);
		for (int i = 1; i < 40; i++)
		{
			radius = this->plys[0].sampleR();
			theta = rand01();
			while (conflict(j, radius, theta))
			{
				//std::cout << "repeat\n";
				radius = this->plys[0].sampleR();
				theta = 2.0f * pi * rand01();
			}
			//std::cout << i << " " << radius << " " << theta << std::endl;
			std::pair<float, float> sample(radius, theta);
			this->plys[j].cross_section_samples.push_back(sample);
		}
	}
}

void Fiber::Yarn::createCrossSectionSamplesUniform()
{
	for (int i = 0; i < 60; i++)
	{
		float radius = this->plys[0].sampleR();
		//while (radius < 0.5f)
		{
			radius = this->plys[0].sampleR();
		}
		//radius = 1.0;
		//float theta = (i / 20.0);
		float theta = (i / 20.0f) - floor(i / 20.0f) + rand01() / 20.0f;
		//float theta = (i / 20.0) - floor(i / 20.0); // rand01();
		//float theta = rand01();
		std::pair<float, float> p(radius, theta);
		crossSectionSamples.push_back(p);
	}
}

void Fiber::Yarn::createCrossSectionSamples()
{
	createPoissonSamples();

	const int ply_num = this->plys.size();
	for (int j = 0; j < ply_num; j++) {
		int i = ply_num - 1 - j;
		std::sort(this->plys[i].cross_section_samples.begin(), this->plys[i].cross_section_samples.end());
		for (int t = 0; t < 20; t++)
			//crossSectionSamples.push_back(std::pair<float, float>(10000.0, 0.0));
			crossSectionSamples.push_back(this->plys[i].cross_section_samples[this->plys[i].cross_section_samples.size() - 1 - t]);
		//for (int t = 0; t < 0; t++)
		//{
		//	crossSectionSamples.push_back(std::pair<float, float>(10000.0, 0.0));
		//}
	}
}


float randomFloat(float min, float max) {
	static std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution<float> distribution(min, max);
	return distribution(generator);
}

ks::vec3 perturbDirection(const ks::vec3& direction, float maxAngleRadians) {
	ks::vec2 uv;
	uv.x() = randomFloat(0.0f, 1.0f);
	uv.y() = randomFloat(0.0f, 1.0f);

	ks::vec3 normalizedDirection = direction.normalized();

	float theta = uv.x() * maxAngleRadians;
	float phi = uv.y() * 2.0f * ks::pi;

	ks::vec3 randomOffset(
		std::sin(theta) * std::cos(phi),
		std::sin(theta) * std::sin(phi),
		std::cos(theta)
	);

	ks::Frame frame(normalizedDirection);

	ks::vec3 perturbedDirection = frame.to_world(randomOffset).normalized();
	return perturbedDirection;
}

float gradientFunction(float t, float r0, float re, float T, float beta, float k, float omega) {
	float A = 1 - std::exp(-beta);
	float gamma = beta * std::exp(-beta) / (1 - std::exp(-beta));
	if (t <= T) {
		return r0 + re * (1 - exp(-beta * t / T));
	}
	else {
		double descending = A * std::exp(-gamma * (t - T) / T) * (1 + k * std::sin(omega * (t - T)));
		return r0 + re * descending;
	}
}

float calculateY(float x, float z0, float t, float ze, float gamma, float alpha = 1.0f) {
	if (x <= t) {
		// Left side of the peak (rising side)
		float scale = x / t;
		return z0 + ze * (1.0f - std::pow(1.0f - scale, alpha));
	}
	else {
		// Right side of the peak (falling side)
		float scale = (x - t) / (1.0f - t);
		return z0 + ze * (1.0f - std::pow(scale, alpha) * gamma);
	}
}

void Fiber::Yarn::simulate_fly_away()
{// Step1: Obtain center of yarn starting point
	const vec3 base_center(0, 0, 0);
	const float base_radius = yarn_radius;
	const int ply_num = this->plys.size();
	// Step2: Sample initial locations of ply-centers in normal plane around starting point
	for (int i = 0; i < ply_num; i++) {
		float angle = 2 * pi * i / ply_num + 0.5 * pi;
		this->plys[i].base_theta = angle;
		this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
	}

	// Step3: Sample initial fiber locations in normal plane around ply-centers using rejection sampling according
	// to the distribution in Sec 4.1
	// RNG rng(rng_seed);
	std::string filename = "data_yarn/fiber_random.txt";
	std::ifstream file(filename);
	std::vector<vec3> randoms(0);
	std::vector<vec2i> pos(0);
	if (!file.is_open()) {
		std::cerr << "Error: Could not open file " << filename << std::endl;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::array<float, 3> floats;
		std::array<int, 2> ints;

		if (ss >> ints[0] >> ints[1] >> floats[0] >> floats[1] >> floats[2]) {
			randoms.push_back(vec3(floats[0], floats[1], floats[2]));
			pos.push_back(vec2i(ints[0], ints[1]));
		}
		else {
			std::cerr << "Error: Incorrect data format in line: " << line << std::endl;
		}
	}

	file.close();

	const int num_of_cores = omp_get_num_procs();

	for (int i = 0; i < ply_num; i++) {
		const int fiber_num = this->plys[i].fibers.size();
#pragma omp parallel for num_threads(num_of_cores) 
		for (int f = 0; f < fiber_num; f++) {
			Fiber& fiber = this->plys[i].fibers[f];
			float radius = this->plys[i].sampleR();
			float theta = 2 * pi * rand01();
			float migration_theta = 2 * pi * rand01();
			fiber.init_radius = radius;
			fiber.init_theta = theta;
			fiber.init_migration_theta = migration_theta;
			fiber.init_vertex = this->plys[i].base_center +
				vec3(radius * std::cosf(theta), radius * std::sinf(theta), 0);
		}

	}

	for (int i = 0; i < pos.size(); i++) {
		vec2i pos_cur = pos[i];
		vec3 random_cur = randoms[i];
		Fiber& fiber = this->plys[pos_cur.x()].fibers[pos_cur.y()];
		if ((pos_cur.x() < ply_num) && (pos_cur.y() < this->plys[0].fibers.size())) {

			fiber.init_radius = random_cur.x();
			fiber.init_theta = random_cur.y();
			fiber.init_migration_theta = random_cur.z();
			fiber.init_vertex =
				this->plys[pos_cur.x()].base_center + vec3(fiber.init_radius * std::cosf(fiber.init_theta),
					fiber.init_radius * std::sinf(fiber.init_theta), 0);
		}
	}

	// Step4: Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center
	// and rotating fiber positions around ply centers

	std::vector<std::vector<std::vector<float>>> rVals(ply_num);
	std::normal_distribution<float> distrbn;
	std::uniform_real_distribution<float> distrbr;

	std::random_device rd;

	std::mt19937 engine(rd());

	// std::cout << this->z_step_num << std::endl;
#pragma omp parallel for num_threads(num_of_cores) 
	for (int i = 0; i < ply_num; i++) {
		const int fiber_num = (int)this->plys[i].fibers.size();
		const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);
		rVals[i].resize(fiber_num);
		for (int f = 0; f < fiber_num; f++) {
			Fiber& fiber = this->plys[i].fibers[f];
			fiber.clear();
			rVals[i][f].clear();
			float prev_mt_local = 0;
			float prev_mt_local_period = 0;

			float loop_radius = -1;

			float beta_loop = plys[i].flyaway_loop_density * this->plys[i].alpha / (this->plys[i].s_i * fiber_num);
			// float beta_loop = 0.165855616;

			// printf("%f\n", beta_loop);

			// roll a dice to decide whether to start a loop period
			if (this->use_flyaways && distrbr(engine) < beta_loop) {
				for (;;) {
					// this sigma is actually twice the sigma of the gaussian distribution
					loop_radius =
						(plys[i].flyaway_loop_r1_mu + 0.5f * plys[i].flyaway_loop_r1_sigma * distrbn(engine)) /
						balance_radius;
					if (loop_radius > (this->plys[i].rho_max * fiber.init_radius))
						break; // to match ours, probably don't clamp it
				}
			}

			for (int step_id = 0; step_id < this->z_step_num; step_id++) {
				const float z = this->z_step_size * (step_id - this->z_step_num / 2.f);
				const float fiber_theta =
					this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
				const float yarn_theta =
					this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
				float local_x, local_y, world_x, world_y;

				// Step5: Vary the distance of cross-sectional fiber positions to their ply center according to fiber
				// migration Sec 4.2

				float migration_theta = fiber.init_migration_theta + this->plys[i].s_i * (fiber_theta);
				float mt_local_period = floor(migration_theta / (2 * pi));
				float mt_local = migration_theta - 2 * pi * mt_local_period;

				if (step_id > 0 && mt_local >= pi && prev_mt_local <= pi && mt_local_period == prev_mt_local_period) {
					// roll a dice to decide whether to start a loop period
					if (this->use_flyaways && distrbr(engine) < beta_loop) {
						for (;;) {
							// this sigma is actually twice the sigma of the gaussian distribution
							loop_radius =
								(plys[i].flyaway_loop_r1_mu + 0.5f * plys[i].flyaway_loop_r1_sigma * distrbn(engine)) /
								balance_radius;
							if (loop_radius > (this->plys[i].rho_max * fiber.init_radius))
								break; // to match ours, probably don't clamp it
						}
					}
					else
						loop_radius = -1;
				}
				prev_mt_local = mt_local;
				prev_mt_local_period = mt_local_period;

				//this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, -1, local_x, local_y);
				this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, loop_radius, local_x, local_y);

				vec3 short_axis = this->plys[i].base_center.normalized(),
					long_axis = vec3(-short_axis.y(), short_axis.x(), 0);
				vec3 local_p = vec3(local_x, local_y, 0.f);
				float _local_x = local_p.dot(short_axis), _local_y = local_p.dot(long_axis);
				_local_x *= this->plys[i].ellipse_short;
				_local_y *= this->plys[i].ellipse_long;

				local_p = _local_x * short_axis + _local_y * long_axis;
				local_x = local_p.x();
				local_y = local_p.y();

				rVals[i][f].push_back(std::sqrt(local_x * local_x + local_y * local_y));

				float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x();
				float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y();

				world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) -
					world_y_before_ply_rotation * std::sinf(yarn_theta);
				world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) +
					world_x_before_ply_rotation * std::sinf(yarn_theta);

				vec3 verIn(world_x, world_y, z);

				//std::cout << verIn.x() << " " << verIn.y() << " " << verIn.z() << std::endl;
				fiber.vertices.push_back(verIn);
			}
		}
	}

	float num_hair_fibers = 0;
	if (this->use_flyaways) {
		std::uniform_real_distribution<float> distrb1;
		std::normal_distribution<float> distrb2;

		const float sig_scale_hair = 0.5f, sig_scale_loop = 0.5f;
		const int min_loop_span = 10;

		float zextent = this->aabb_micro_ct.pMax.z() - this->aabb_micro_ct.pMin.z();
		const float base_radius = this->yarn_radius;
		for (int i = 0; i < this->getPlyNum(); ++i) {

			float angle = 2 * pi * i / this->getPlyNum();
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
			const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);

			int nloop = static_cast<int>(std::floor(plys[i].flyaway_loop_density * zextent + 0.5f));
			int loop_hair = 0;
			if (0 && (nloop > 0)) {
				std::vector<vec2i> locs;
				int fiber_num = static_cast<int>(this->plys[i].fibers.size());
				for (int j = 0; j < fiber_num; ++j) {
					const Fiber& curFiber = this->plys[i].fibers[j];
					int totVtx = static_cast<int>(curFiber.vertices.size());
					for (int k = 1; k + 1 < totVtx; ++k)
						if (rVals[i][j][k] > rVals[i][j][k - 1] && rVals[i][j][k] > rVals[i][j][k + 1]) {
							locs.push_back(vec2i(j, k));
						}
				}

				std::shuffle(locs.begin(), locs.end(), engine);

				std::vector<std::vector<vec3>> loop_fibers;
				for (int j = 0; j < nloop && j < static_cast<int>(locs.size()); ++j) {
					int fid = locs[j][0];
					const std::vector<float> curRs = rVals[i][fid];
					Fiber& curFiber = this->plys[i].fibers[fid];
					int totVtx = static_cast<int>(curFiber.vertices.size());
					const int k = locs[j][1];
					int k0 = std::max(k - min_loop_span, 0);
					int k1 = std::min(k + min_loop_span, totVtx - 1);
					while (k0 > 0 && curRs[k0 - 1] < curRs[k0])
						--k0;
					while (k1 + 1 < totVtx && curRs[k1 + 1] < curRs[k1])
						++k1;

					float r1;
					for (;;) {
						// this sigma is actually twice the sigma of the gaussian distribution
						r1 = plys[i].flyaway_loop_r1_mu +
							sig_scale_loop * plys[i].flyaway_loop_r1_sigma * distrb2(engine);
						if (r1 > 1.f * curRs[k])
							break; // to match ours, probably don't clamp it
					}

					float ratio = r1 / curRs[k];
					Fiber fiber;
					fiber.vertices.resize(0);
					for (int t = k0 + 1; t <= k; ++t) {
						// float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k0)/static_cast<float>(k - k0);
						float v = 1.0f + (ratio - 1.0f) * std::sin(0.5f * pi * static_cast<float>(t - k0) /
							static_cast<float>(k - k0));
						vec3 point = curFiber.vertices[t];
						point.x() *= v;
						point.y() *= v;
						fiber.vertices.push_back(point);

					}
					for (int t = k1 - 1; t > k; --t) {
						// float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k1)/static_cast<float>(k - k1);
						float v = 1.0f + (ratio - 1.0f) * std::sin(0.5f * pi * static_cast<float>(t - k1) /
							static_cast<float>(k - k1));
						vec3 point = curFiber.vertices[t];
						point.x() *= v;
						point.y() *= v;
						fiber.vertices.push_back(point);
					}
					if (fiber.vertices.size() > 3) {
						loop_hair++;
						this->plys[i].fibers.push_back(fiber);
					}
				}
			}

			plys[i].num_hair_fibers = 0;

			int nhair = static_cast<int>(std::floor(plys[i].flyaway_hair_density * zextent + 0.5f));
			for (int j = 0; j < nhair;) {
				float hairtype = randomFloat(0, 1);
				Fiber fiber;
				float z0 = this->aabb_micro_ct.pMin.z() + distrb1(engine) * zextent;

				float ze;
				for (;;) {
					ze = plys[i].flyaway_hair_ze_mu +
						sig_scale_hair * plys[i].flyaway_hair_ze_sigma *
						distrb2(engine); //+ plys[i].flyaway_hair_ze_sigma*distrb1(engine);
					if (ze > 0)
						break;
				}

				//					float ze = plys[i].flyaway_hair_ze_mu +
				//plys[i].flyaway_hair_ze_sigma*distrb1(engine);

				float r0;
				for (;;) {
					r0 = plys[i].flyaway_hair_r0_mu + sig_scale_hair * plys[i].flyaway_hair_r0_sigma * distrb2(engine);
					if (r0 >= 0) {
						r0 = std::min(r0, balance_radius);
						break;
					}
				}

				float re = 0;
				for (;;) {
					re = plys[i].flyaway_hair_re_mu + sig_scale_hair * plys[i].flyaway_hair_re_sigma * distrb2(engine);
					if (re > 0)
						break;
				}

				float p0 = 2.0f * pi * distrb1(engine);
				float pe =
					plys[i].flyaway_hair_pe_mu + sig_scale_hair * plys[i].flyaway_hair_pe_sigma * distrb2(engine);

				int nstep = 24;
				std::vector<vec3> vars(0);

				if (hairtype > 0.8) {
					// parabolic type

					float t = randomFloat(0.7, 0.95);   // x-coordinate of the peak
					float alpha = randomFloat(1.4, 3); // Asymmetry factor
					float gamma = randomFloat(0.0, 0.2); // descent factor

					for (int k = 0; k <= nstep; ++k) {
						vec3 cur;
						// cur[0] = r0_e + re_e * static_cast<float>(k) / static_cast<float>(nstep);
						// cur[1] = z0_e + ze_e * static_cast<float>(k) / static_cast<float>(nstep);
						// cur[2] = p0_e + pe_e * static_cast<float>(k) / static_cast<float>(nstep);

						float x = static_cast<float>(k) / static_cast<float>(nstep);
						cur[0] = calculateY(x, r0, t, re, gamma, alpha);
						cur[1] = z0 + ze * static_cast<float>(k) / static_cast<float>(nstep);

						// p0 + cur[1] * 2 * Zhao::pi / this->plys[i].alpha; //
						cur[2] = p0 + pe * static_cast<float>(k) / static_cast<float>(nstep);

						//float diff_ze = plys[i].flyaway_loop_r1_sigma * distrb2(engine); //+ plys[i].flyaway_hair_ze_sigma*distrb1(engine);
						//ze += diff_ze;

						//float diff_re = 0;
						//for (;;) {
						//	diff_re = plys[i].flyaway_loop_r1_sigma * distrb2(engine);
						//	if (diff_re + re > 0) {
						//		re = diff_re + re;
						//		break;
						//	}
						//}

						//float diff_pe = plys[i].flyaway_loop_r1_sigma * distrb2(engine);
						//pe += diff_pe;

						vars.push_back(cur);
					}

					//std::cout << " r0 is " << r0 << " re is " << re << std::endl;
					/* Creating fiber curve */

					bool visible = false;
					for (int k = 0; k <= nstep; ++k) {
						const vec3& cur = vars[k];
						vec3 pos;
						pos[0] = cur[0] * std::cos(cur[2]);
						pos[1] = cur[0] * std::sin(cur[2]);
						pos[2] = cur[1];
						// Crop flyaway fibers using the ply's bounding box
						if ((pos[2] < this->aabb_micro_ct.pMin.z()) || (pos[2] > this->aabb_micro_ct.pMax.z()))
							break;

						if ((cur[0] < balance_radius) && (!visible)) {
							visible = true;
							fiber.vertices.resize(0);
						}
						fiber.vertices.push_back(pos + this->plys[i].base_center);
					}
					if (fiber.vertices.size() > 3) {
						plys[i].num_hair_fibers++;
						this->plys[i].fibers.push_back(fiber);
						++j;
					}
				}
				else if (hairtype > 0.5) {

					// ALL Random direction start from first step

					vec3 p_cur;
					p_cur[0] = r0 * std::cos(p0);
					p_cur[1] = r0 * std::sin(p0);
					p_cur[2] = z0;
					vec3 dir;
					float step_size = 0;
					vars.push_back(p_cur);
					{
						vec3 cur;
						cur[0] = r0 + re;
						cur[1] = z0 + ze;

						// p0 + cur[1] * 2 * Zhao::pi / this->plys[i].alpha; //
						cur[2] = p0 + pe;

						vec3 p0;
						p0[0] = cur[0] * std::cos(cur[2]);
						p0[1] = cur[0] * std::sin(cur[2]);
						p0[2] = cur[1];
						dir = (p0 - p_cur);
						step_size = dir.norm() / nstep;
						dir = dir.normalized();
					}

					for (int k = 1; k <= nstep; k++) {
						dir = perturbDirection(dir, 0.5);
						p_cur = p_cur + step_size * dir;
						vars.push_back(p_cur);
					}

					for (int k = 0; k <= nstep; ++k) {
						const vec3& cur = vars[k];
						vec3 pos = cur;
						// Crop flyaway fibers using the ply's bounding box
						if ((pos[2] < this->aabb_micro_ct.pMin.z()) || (pos[2] > this->aabb_micro_ct.pMax.z()))
							break;

						fiber.vertices.push_back(pos + this->plys[i].base_center);
					}
					if (fiber.vertices.size() > 3) {
						plys[i].num_hair_fibers++;
						this->plys[i].fibers.push_back(fiber);
						++j;
					}
				}
				else {
					// Shuang's version

					for (int k = 0; k <= nstep; ++k) {
						vec3 cur;
						// cur[0] = r0_e + re_e * static_cast<float>(k) / static_cast<float>(nstep);
						// cur[1] = z0_e + ze_e * static_cast<float>(k) / static_cast<float>(nstep);
						// cur[2] = p0_e + pe_e * static_cast<float>(k) / static_cast<float>(nstep);

						cur[0] = r0 + re * static_cast<float>(k) / static_cast<float>(nstep);
						cur[1] = z0 + ze * static_cast<float>(k) / static_cast<float>(nstep);

						// p0 + cur[1] * 2 * Zhao::pi / this->plys[i].alpha; //
						cur[2] = p0 + pe * static_cast<float>(k) / static_cast<float>(nstep);

						//float diff_ze = plys[i].flyaway_loop_r1_sigma * distrb2(engine); //+ plys[i].flyaway_hair_ze_sigma*distrb1(engine);
						//ze += diff_ze;

						//float diff_re = 0;
						//for (;;) {
						//	diff_re = plys[i].flyaway_loop_r1_sigma * distrb2(engine);
						//	if (diff_re + re > 0) {
						//		re = diff_re + re;
						//		break;
						//	}
						//}

						//float diff_pe = plys[i].flyaway_loop_r1_sigma * distrb2(engine);
						//pe += diff_pe;

						vars.push_back(cur);
					}

					/* Creating fiber curve */

					bool visible = false;
					for (int k = 0; k <= nstep; ++k) {
						const vec3& cur = vars[k];
						vec3 pos;
						pos[0] = cur[0] * std::cos(cur[2]);
						pos[1] = cur[0] * std::sin(cur[2]);
						pos[2] = cur[1];
						// Crop flyaway fibers using the ply's bounding box
						if ((pos[2] < this->aabb_micro_ct.pMin.z()) || (pos[2] > this->aabb_micro_ct.pMax.z()))
							break;

						if ((cur[0] < balance_radius) && (!visible)) {
							visible = true;
							fiber.vertices.resize(0);
						}
						fiber.vertices.push_back(pos + this->plys[i].base_center);
					}
					if (fiber.vertices.size() > 3) {
						plys[i].num_hair_fibers++;
						this->plys[i].fibers.push_back(fiber);
						++j;
					}
				}
			}
			plys[i].num_hair_fibers += loop_hair;
			num_hair_fibers += plys[i].num_hair_fibers;
		}
		std::cout << "we have " << num_hair_fibers << " flyaway" << std::endl;
	}
	std::cout << "Checking..." << std::endl;
	std::cout << num_hair_fibers << std::endl;
	int bad_count = 0;
	for (int i = 0; i < this->plys.size(); i++) {
		const int fiberNum = (int)this->plys[i].fibers.size();
		std::cout << "ply " << i << " is " << fiberNum << std::endl;
		for (int f = 0; f < fiberNum; f++) {
			const int vertexNum = (int)this->plys[i].fibers[f].vertices.size();
			// std::cout << i << " " << f << " has " << vertexNum << std::endl;
			for (int v = 1; v < vertexNum - 1; v++) {
				vec3 prev = this->plys[i].fibers[f].vertices[v - 1];
				vec3 curr = this->plys[i].fibers[f].vertices[v];
				vec3 next = this->plys[i].fibers[f].vertices[v + 1];
				vec3 dir1 = (next - curr).normalized();
				vec3 dir2 = (curr - prev).normalized();
				if (dir1.dot(dir2) < 0.5) {
					bad_count++;
				}
			}
		}
	}
	if (bad_count > 0)
		std::cout << "Bad count = " << bad_count << std::endl;
	//std::ofstream outFile("D:/output_RT.txt");

	//if (!outFile.is_open()) {
	//	std::cerr << "bad!" << std::endl;
	//}

	//for (int i = 0; i < this->plys.size(); i++) {
	//	const int fiberNum = (int)this->plys[i].fibers.size();
	//	for (int f = 0; f < fiberNum; f++) {
	//		const int vertexNum = (int)this->plys[i].fibers[f].vertices.size();
	//		// std::cout << i << " " << f << " has " << vertexNum << std::endl;
	//		for (int v = 0; v < vertexNum; v++) {
	//			outFile << this->plys[i].fibers[f].vertices[v].x() << " "
	//				<< this->plys[i].fibers[f].vertices[v].y() << " " << this->plys[i].fibers[f].vertices[v].z() << std::endl;
	//		}
	//	}
	//}
	//outFile.close();
}