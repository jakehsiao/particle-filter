/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;	

  	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i<num_particles; i++){
	  Particle this_particle;
	  this_particle.x = dist_x(gen);
	  this_particle.y = dist_y(gen);
	  this_particle.theta = dist_theta(gen);
	  this_particle.id = i;
	  this_particle.weight = 1;
	  weights.push_back(1);
	  particles.push_back(this_particle);
	  
	}
	is_initialized = true;
	

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0, std_pos[0]);
	std::normal_distribution<double> dist_y(0, std_pos[1]);
	std::normal_distribution<double> dist_theta(0, std_pos[2]);

	int pc = 0;
	for (int i = 0; i < num_particles; i++){
	  
	  // add noise
	  Particle* this_particle = &particles[i];
	  this_particle->x += dist_x(gen);
	  this_particle->y += dist_y(gen);
	  this_particle->theta += dist_theta(gen);

	  // add motion
	  float theta = this_particle->theta;
	  float new_theta = theta + yaw_rate * delta_t;
	  float v_over_yr = velocity / yaw_rate;
	  this_particle->x += v_over_yr * (sin(new_theta) - sin(theta));
	  this_particle->y += v_over_yr * (cos(theta) - cos(new_theta));
	  this_particle->theta = new_theta;
	  	
	} 

}

double** ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	int m = predicted.size();
	int n = observations.size();
	double** errors = new double*[m];

	for (int i = 0; i < m; i++){
	  errors[i] = new double[2];

	  LandmarkObs lm = predicted[i];
	  //std::cout << "pred:" << lm.x << " " << lm.y;
	  double best_distance = 9999;
	  LandmarkObs obs_bestfit;

	  // get the nearest neighbour
	  for (int j=0; j < n; j++){
	    LandmarkObs obs = observations[j];
	    double delta_x = lm.x - obs.x;
	    double delta_y = lm.y - obs.y;
	    double distance = sqrt(delta_x * delta_x + delta_y * delta_y);
	    //update the nearest neibour
	    if (distance < best_distance){
	      obs_bestfit = obs;
	      best_distance = distance;
	    }
	  }
	  //std::cout << " obs:" << obs_bestfit.x << " " << obs_bestfit.y;

	  errors[i][0] = obs_bestfit.x - lm.x;
	  errors[i][1] = obs_bestfit.y - lm.y;

	  //std::cout << " err:" << dist(obs_bestfit.x, obs_bestfit.y, lm.x, lm.y) << std::endl;
	}

	return errors;



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
  	int n = observations.size(); // num of observations

  	for (int i=0; i<num_particles; i++){
	  Particle this_particle = particles[i];
	  std::vector<LandmarkObs> predicted;

	  for (int j=0; j<map_landmarks.landmark_list.size(); j++){
	    Map::single_landmark_s this_lm_map = map_landmarks.landmark_list[j];
	    LandmarkObs this_lm;
	    this_lm.x = this_lm_map.x_f;
	    this_lm.y = this_lm_map.y_f;

	    double dist_from_car = dist(this_particle.x, this_particle.y, this_lm.x, this_lm.y);
	    if (dist_from_car <= sensor_range){
	      predicted.push_back(this_lm);
	    }
	  }
	  


	  std::vector<LandmarkObs> observations_transformed;
	  for (int j=0; j<n; j++){
	    LandmarkObs this_obs = observations[j];
	    LandmarkObs transformed;

	    double cos_theta = cos(this_particle.theta);
	    double sin_theta = sin(this_particle.theta);
	    transformed.x = cos_theta * this_obs.x + sin_theta * this_obs.y;
	    transformed.y = sin_theta * this_obs.x - cos_theta * this_obs.y;

	    transformed.x += this_particle.x;
	    transformed.y += this_particle.y;

	    observations_transformed.push_back(transformed);


	  }
	  //std::cout << std::endl;
	  
	  //  compare to the observation of it, returns the errors
	  double** errors = dataAssociation(predicted, observations_transformed);
	  //  update the weight according to errors
	  int m = predicted.size();
	  double* w = &weights[i];
	  *w = 1;
	  for (int j=0; j<m; j++){
	    LandmarkObs this_p_lm = predicted[j];
	    double distance_from_car = dist(this_particle.x, this_particle.y, this_p_lm.x, this_p_lm.y);
	    double params_std = 1/distance_from_car; // TOTUNE
	    double std_x = std_landmark[0] * distance_from_car * params_std;
	    double std_y = std_landmark[1] * distance_from_car * params_std;
	    double var_x = std_x * std_x; 
	    double var_y = std_y * std_y;
	    double sqr_delta_x = pow(errors[j][0], 2);
	    double sqr_delta_y = pow(errors[j][1], 2);

	    // weight calc 1: by gaussian prob density
	    *w *= exp((-0.5) * (sqr_delta_x / var_x + sqr_delta_y / var_y)) / sqrt(2 * 3.14159 * var_x * var_y);
	    // weight calc 2: by eu dist
	    //*w += sqrt(sqr_delta_x + sqr_delta_y) / m;

	  }
	  particles[i].weight = *w;
	//std::cout << "weight: " << *w << std::endl << std::endl;
	}
	//std::cout << "FINISH ONE ITER" << std::endl;

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine generator;
	std::discrete_distribution<int> distribution (weights.begin(), weights.end());
	std::vector<Particle> old_particles = particles;
	// std::vector<Particle> *new_particles = new std::vector<Particle>;

	int new_index;
	for (int i = 0; i < num_particles; i++){
	  new_index = distribution(generator);
	  particles[i] = old_particles[new_index];
	}

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
