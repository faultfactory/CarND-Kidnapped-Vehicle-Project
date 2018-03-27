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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles=50;
	default_random_engine ranGen;
	particles.resize(num_particles);
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);
	int	p_id= 0; 
	for(std::vector<Particle>::iterator iter_p = particles.begin(); iter_p != particles.end(); iter_p++){
		//set particle Id. 
		iter_p->id=p_id; 
		iter_p->x=dist_x(ranGen);
		iter_p->y=dist_y(ranGen);
		iter_p->theta=dist_theta(ranGen);
		iter_p->weight=1.0;
	}
	is_initialized=true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine ranGen;


	for(std::vector<Particle>::iterator iter_p = particles.begin(); iter_p != particles.end() ; iter_p++){
		double tmp_x, tmp_y, tmp_theta; 
		double px = iter_p->x;
		double py = iter_p->y;
		double theta = iter_p->theta; 
		// Borrowing some code from the UKF project to protect against YawRate==0
		// discarding the yaw rate accel
		
		if(fabs(yaw_rate)>0.0001){
			tmp_x=px+(velocity/yaw_rate)*(sin(theta+yaw_rate*delta_t)-sin(theta));
			tmp_y=py+(velocity/yaw_rate)*(-cos(theta+yaw_rate*delta_t)+cos(theta));
		}
		else{
			tmp_x=px+(velocity*cos(theta)*delta_t);
			tmp_y=py+(velocity*sin(theta)*delta_t);
			cout<<"yawrateLow"<<endl;
		}
		tmp_theta= theta + yaw_rate*delta_t;
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
	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGrathetacs/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	//Convert the landmarks into an observation vector because we must submit the .cpp file 
	//and can't modify the header according to documentation. (C'mon guys)
	std::vector<LandmarkObs> landmarks;
	for(auto iter_lm = map_landmarks.landmark_list.begin(); iter_lm!=map_landmarks.landmark_list.end();iter_lm++){
		LandmarkObs lm;
		lm.x=iter_lm->x_f; 
		lm.y=iter_lm->y_f; 
		lm.id=iter_lm->id_i; 
		landmarks.push_back(lm);
	}
	
	// Compute some constants for use later:
	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	// Gassuian norm used in weight calculation; 
	double gs_nrm = (1.0/(2.0*M_PI*sig_x*sig_y));
	// Loop through each particle
	for(auto iter_p = particles.begin(); iter_p!=particles.end();iter_p++){
		//Clear storage values for this particle; 
		iter_p->associations.clear();
		iter_p->sense_x.clear();
		iter_p->sense_y.clear();

		//Create a vector of new observations;
		std::vector<LandmarkObs> pred_obs;
		// For each of the observations;
		for(std::vector<LandmarkObs>::const_iterator iter_obs = observations.begin(); iter_obs != observations.end();iter_obs++){
			//create a new obs
			LandmarkObs obs_out;
			//create rotated observation values
			double x_rot = (iter_obs->x * cos(iter_p->theta)) + (iter_obs->y * sin(iter_p->theta));
			double y_rot = -(iter_obs->x * sin(iter_p->theta)) + (iter_obs->y * cos(iter_p->theta));
			//add rotated observation to particle position and assign to temporary observation 
			obs_out.x = iter_obs->x+x_rot; 
			obs_out.y = iter_obs->y+y_rot; 
			// push to output vector
			pred_obs.push_back(obs_out);
		}
		
		// Execute the data association. The structure of that helper doesn't work for me so i'm ignoring it. 
		// For each landmark 
		for(auto iter_pobs = pred_obs.begin(); iter_pobs !=pred_obs.end();iter_pobs++){
			// initialize score variable to sensor range as huge starting value; 
			double dist_score = sensor_range;
			int lm_asc=0;

			for(auto iter_lm = landmarks.begin();iter_lm != landmarks.end(); iter_lm++){
				// calculate distance to the landmark to reject detections beyond sensor range
				double lm_dist=dist(iter_p->x,iter_p->y,iter_lm->x,iter_lm->y);
				if(lm_dist>sensor_range){
					//Do nothing
				}
				else{
					double od=dist(iter_pobs->x,iter_pobs->y,iter_lm->x,iter_lm->y);
					// if the observation distance is smaller than the distance currently recorded
					if (dist_score>od){
						// save the new distance
						dist_score=od;
						// and the landmark id. 
						lm_asc = iter_lm->id;
					}
				}
			}
			if(lm_asc !=0){
				iter_p->associations.push_back(lm_asc);
				iter_p->sense_x.push_back(iter_pobs->x);
				iter_p->sense_y.push_back(iter_pobs->y);
			}
		}
	}
	weights.clear();
	// start a new particle loop for clarity. 
	auto lm_list=map_landmarks.landmark_list;
	for(auto iter_p = particles.begin(); iter_p!=particles.end();iter_p++){
		// check if any data was associated with this particle.
		if(iter_p->associations.empty()){
			weights.push_back(0.0);
		} 
		else {
			// new storage for particle weight;
			double wt=1.0;
			//grab vectors with data
			auto asc=iter_p->associations;
			auto sx=iter_p->sense_x;
			auto sy=iter_p->sense_y;
			cout<<"number of asc "<<asc.size()<<endl;
			for(int i = 0; i<asc.size();i++){
				// sub 1 from the asc index to get the correct value; 
				auto lm=lm_list[asc[i]-1];
				double u_x = lm.x_f;
				double u_y = lm.y_f;
				double x = sx[i];
				double y = sy[i];		
				double expn = ((x-u_x)*(x-u_x)/(2.0*sig_x*sig_x))+((y-u_y)*(y-u_y)/(2.0*sig_y*sig_y));
				wt=wt*gs_nrm*exp(-expn);
				cout<<"WT: "<<wt<<endl;
			}
			iter_p->weight=wt;
			weights.push_back(wt);

		}
	}
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

		std::vector<Particle> new_p; 
	default_random_engine gen(rand());
	std::discrete_distribution<> dist_p(std::begin(weights), std::end(weights));
	double chsn_indx = dist_p(gen);
	for(int i=0; i < num_particles;i++){
	new_p.push_back(particles[chsn_indx]);
	}

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
