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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "map.h"

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine ranGen;
	num_particles=1000;	
	std::vector<Particle> particles(1000);
	
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);

	int	p_id= 0; 

	// Practicing standard library functions. 
	for(std::vector<Particle>::iterator iter_p = particles.begin(); iter_p != particles.end() ; iter_p ++){
		//set particle Id. 
		iter_p->id=p_id; 
		p_id++;
		iter_p->x=dist_x(ranGen);
		iter_p->y=dist_y(ranGen);
		iter_p->theta=dist_theta(ranGen);
		iter_p->weight=1.0;
		}


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//create temporary pose variables
	double tmp_x; double tmp_y; double tmp_theta;

	default_random_engine ranGen;
	for(std::vector<Particle>::iterator iter_p = particles.begin(); iter_p != particles.end() ; iter_p ++){
		tmp_x= iter_p->x + (velocity/yaw_rate)*(sin((iter_p->theta + yaw_rate*delta_t) - sin(iter_p->theta)));
		tmp_y= iter_p->y + (velocity/yaw_rate)*(cos((iter_p->theta) - cos(iter_p->theta + yaw_rate*delta_t)));
		tmp_theta= iter_p->theta + yaw_rate;

		// TODO: If time permits consider doing this outside of for loop with  zero center and 
		// Would just involve adding at the end of each equation instead of making a distribution. 
		normal_distribution<double> dist_x(tmp_x,std_pos[0]);
		normal_distribution<double> dist_y(tmp_y,std_pos[1]);
		normal_distribution<double> dist_theta(tmp_theta,std_pos[2]);
		// Replace existing position values with those including travel and Guassian noise added; 
		iter_p->x=dist_x(ranGen);
		iter_p->y=dist_y(ranGen);
		iter_p->theta=dist_theta(ranGen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Cycle through all the available landmark predictions 
	for(std::vector<LandmarkObs>::iterator iter_pred = predicted.begin(); iter_pred !=predicted.end();iter_pred++){
		double dist_score; // initializing score value
		std::vector<LandmarkObs>::iterator bestMatchPtr; // pointer to best match observation
		for(std::vector<LandmarkObs>::iterator iter_obs = observations.begin(); iter_obs != observations.end() ; iter_obs ++){
			//Set values for start
			if(iter_pred==predicted.begin()){
				dist_score = dist(iter_obs->x,iter_obs->y,iter_pred->x,iter_pred->y);
				bestMatchPtr = iter_obs;
			}
			//Compare on every other cycle and overwrite if closer
			else if(dist_score < dist(iter_obs->x,iter_obs->y,iter_pred->x,iter_pred->y)){
					dist_score = dist(iter_obs->x,iter_obs->y,iter_pred->x,iter_pred->y);
					bestMatchPtr = iter_obs;
			}
		}
		//set observation id to the correct landmark 
		bestMatchPtr->id = iter_pred->id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	// convert observations to map coordinates

	for(std::vector<Particle>::iterator iter_p = particles.begin(); iter_p != particles.end() ; iter_p++){
		std::vector<LandmarkObs> p_obs;
		std::vector<LandmarkObs> obs_free = observations; // observations vector is const.. can't modify it. 
		for(std::vector<Map::single_landmark_s>::const_iterator iter_lm = map_landmarks.landmark_list.begin(); iter_lm !=map_landmarks.landmark_list.end();iter_lm++){
			LandmarkObs obs;
			//Generate predicted measurement for each location
			obs.id = iter_lm->id_i;		
			double map_x,map_y; // temp location varaibles for predicted measurement in map coordinates
			map_x = iter_lm->x_f -iter_p->x;
			map_y = iter_lm->y_f -iter_p->y;
			// Execute coordinate tranformation using the particle orientation value. 
			obs.x = map_x*cos(iter_p->theta)+map_y*sin(iter_p->theta);
			obs.y = map_y*cos(iter_p->theta)-map_x*sin(iter_p->theta);
			// if the landmark distance exists within the range of the sensor
			if(dist(0,obs.x,0,obs.y)<sensor_range){
				// Push the predicted observation into the new vector. 
				p_obs.push_back(obs);
			}
		}
		//Run data association and determine which observations identify which landmark
		dataAssociation(p_obs,obs_free);
		// iterate through each landmark in p_obs and determine which of these should not be included 
		double weight = 1; 
		for(std::vector<LandmarkObs>::iterator iter_pbs=p_obs.begin(); iter_pbs != p_obs.end(); iter_pbs++){
			for(std::vector<LandmarkObs>::iterator iter_obs=obs_free.begin(); iter_obs!=obs_free.end();iter_obs++){
			}

		}
	}



}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
