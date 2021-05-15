/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to --> OK
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if (is_initialized) {
    return;
  }
  num_particles = 100;  // TODO: Set the number of particles
  std::default_random_engine generator;
  std::normal_distribution<double> distribution_x(x,std[0]);
  std::normal_distribution<double> distribution_y(y,std[1]);
  std::normal_distribution<double> distribution_theta(theta,std[2]);

  for (int i = 0; i < num_particles; ++i) {
  Particle particle;
  particle.id = i;
  particle.x = distribution_x(generator);
  particle.y = distribution_y(generator);
  particle.theta = distribution_theta(generator);

  particle.weight = 1.0;
    particles.push_back(particle);
  }
  is_initialized = true;

}

void ParticleFilter::prediction (double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise. --> OK
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine generator;
  std::normal_distribution<double> distribution_x(0,std_pos[0]);
  std::normal_distribution<double> distribution_y(0,std_pos[1]);
  std::normal_distribution<double> angle_theta(0,std_pos[2]);
  for (int i = 0; i < num_particles; ++i) {
    if (fabs(yaw_rate) <= 0.001) {      
      particles[i].x = particles[i].x + velocity * delta_t * cos(particles [i].theta);
      particles[i].y = particles[i].x + velocity * delta_t * sin(particles [i].theta);   

  }
    else {
          particles [i].x = particles [i].x + ( velocity / yaw_rate)* (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles [i].theta));
    particles [i].y = particles [i].y + ( velocity / yaw_rate)* (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
    particles [i].theta = particles [i].theta + yaw_rate * delta_t;
                    
    } 
  
    particles [i].x = particles [i].x + distribution_x(generator);
    particles [i].y = particles [i].y + distribution_y(generator);
    particles [i].theta = particles [i].theta + angle_theta(generator);

  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
int m_observation = observations.size();
int m_prediction = predicted.size();
  
    for(int i=0; i<m_observation; i++){
        double observations_x= observations[i].x;
        double observations_y= observations[i].y;
        int min_id=-1;
        double min_dist=9999999.9;

    for(int j=0; j<m_prediction; j++) {
        double predicted_x= predicted[j].x;
        double predicted_y= predicted[j].y;
        double distance= dist(predicted_x, predicted_y, observations_x, observations_y);
         if(distance<min_dist){
                min_dist=distance;
                min_id= predicted[j].id;
      }

}
         observations[i].id= min_id;

}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
 double stdlm_x = std_landmark[0];
 double stdlm_y = std_landmark[1];
  for (int i = 0; i < num_particles; i++) {
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    std::vector<LandmarkObs> validLandmarks;
    for (int j=0; j < map_landmarks.landmark_list.size();j++) {
      float landmarkX = map_landmarks.landmark_list[j].x_f;
      float landmarkY = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dX = x - landmarkX;
      double dY = y - landmarkY;
      double distance = dist(x, y, dX, dY);
      if (distance<=sensor_range){
        validLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY});
      
    }
  }
 std::vector<LandmarkObs>mappedObservations;
  for (int j = 0; j< observations.size(); j++) {
    double x = particles[j].x;
    double y = particles[j].y;
    double theta = particles[j].theta;
    vector<LandmarkObs> validLandmarks;
    double xm = x + cos(theta) * observations[j].x - sin(theta) * observations[j].y;
    double ym = y + sin(theta) * observations[j].x + cos(theta) * observations[j].y;
    mappedObservations.push_back(LandmarkObs { observations[j].id, xm, ym});
  }
    dataAssociation(validLandmarks, mappedObservations);
   double particles_weight = 1.0;
    
    for (int j=0; j < mappedObservations.size(); j++){
      double observations_x = mappedObservations[j].x;
      double observations_y = mappedObservations[j].y;
      int landmarkID = mappedObservations[j].id;
      double landmarkX, landmarkY;
      
    for (int k=0; k < validLandmarks.size(); k++) 
    {
      if (validLandmarks[k].id == landmarkID) {
        
        landmarkX = validLandmarks[k].x;
        landmarkY = validLandmarks[k].y;
      }
    }
    double dX = observations_x - landmarkX;
    double dY = observations_y - landmarkY;
      
      
    double weight = (1 / (2 * M_PI*stdlm_x*stdlm_y)) * exp(-(dX*dX / (2 * stdlm_x*stdlm_y) + (dY*dY / (2 * stdlm_x*stdlm_y))));
if (weight == 0) {
  particles_weight = particles_weight *0.00001;
}
      else{
        particles_weight = particles_weight * weight;
      }
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
 vector<double> p_weights;
   for (int i = 0; i < particles.size(); i++) {
     p_weights.push_back(particles[i].weight);
   }
    std::vector<Particle> new_particles;
    /* Random number engine class that generates pseudo random numbers */
    std::default_random_engine gen;
    /*std::discrete_distribution produces random integers on the interval [0, n), where the probability 
    of each individual integer i is defined as w i/S, that is the weight of the ith integer divided by the sum of all n weights.*/
    std::discrete_distribution<size_t> distr_index(p_weights.begin(), p_weights.end());
    /* Create new particles with probability proportional to their weight */
    for (auto i = 0; i < particles.size(); i++) {
        new_particles.push_back(particles[distr_index(gen)]);
    }
    /* Copy it to the original particle vector */
    particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}