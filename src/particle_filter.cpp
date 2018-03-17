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
#include <assert.h>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y,    std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for(int i = 0; i < num_particles; i++) {
        Particle *p = new Particle();
        p->id = i;
        p->x = dist_x(gen);
        p->y = dist_y(gen);
        p->theta = dist_theta(gen);
        particles.push_back(*p);

    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0,    std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    for(auto & p : particles) {


        if (fabs(yaw_rate) < 0.00001) {
             p.x += velocity * delta_t * cos(p.theta);
             p.y += velocity * delta_t * sin(p.theta);
           }
       else {
            p.x = p.x + (velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
            p.y = p.y + (velocity/yaw_rate)*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t));
            p.theta = p.theta + yaw_rate * delta_t;
        }


        p.x = p.x + dist_x(gen);
        p.y = p.y + dist_y(gen);
        p.theta = p.theta + dist_theta(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    std::vector<double> d;

    for(auto & p: observations) {
        d.clear();
        std::transform(predicted.begin(),predicted.end(),std::back_inserter(d),
                        [&p] (LandmarkObs  a) {return dist(p.x,p.y,a.x,a.y);});
        auto m = std::min_element(d.begin(),d.end());
        int idx = m - d.begin();
        auto closest = predicted.at(idx);
//        cout << p.x << "," << p.y << "/" << closest.id <<endl;
        p.id = closest.id;
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


    std::vector<LandmarkObs> transformedLOs;
    std::vector<LandmarkObs> allPredictions;
    std::vector<LandmarkObs> predictions;

    // Convenience
    const double s_x = std_landmark[0];
    const double s_y = std_landmark[1];

    // Transform the map objects to landmark objects
    // We do this once and filter for particles later on
    //
    std::transform(map_landmarks.landmark_list.begin(),
                   map_landmarks.landmark_list.end(),
                   std::back_inserter(allPredictions),
                   [](Map::single_landmark_s lm) {
        LandmarkObs no {
        lm.id_i,
        lm.x_f,
        lm.y_f};
        return no;
    });


    // For each particle
    for(auto & p : particles) {
        // Clear particle specific lists
        //
         transformedLOs.clear();
         predictions.clear();
        // Create a list of transformed observations
        //

        std::transform(observations.begin(),
                       observations.end(),
                       std::back_inserter(transformedLOs),

                       [&p](LandmarkObs lo) {
                        LandmarkObs ne{
                        lo.id,
                        lo.x*cos(p.theta) - lo.y*sin(p.theta)+p.x,
                        lo.x*sin(p.theta) + lo.y*cos(p.theta) +p.y};
                        return ne;
        });

        // Get a list of all predictions that are within sensor range.
        // Note that we are using a simple comparison (i.e. a rectangle), and
        // not a precise distance calculation. This will be more efficient in
        // case that we do not have a very high number of landmarks
        //
        std::copy_if(allPredictions.begin(), allPredictions.end(), std::back_inserter(predictions),
                     [&p,&sensor_range] (const LandmarkObs lo) {
                         return (fabs(lo.x-p.x) <= sensor_range) && (fabs(lo.y-p.y) <= sensor_range);
                     });

        dataAssociation(predictions,transformedLOs);

        // Reset weight:
        p.weight = 1.0;
        for(auto const& tlo : transformedLOs) {

            // Find prediction that matches id
            //
            auto associated_prediction =  * std::find_if(predictions.begin(),predictions.end(),
                          [&tlo](LandmarkObs lobby) {
//                        cout << lobby.id << ":" << tlo.id << ":" << (lobby.id == tlo.id) << endl;
                        return lobby.id == tlo.id;
            }) ;


            double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(associated_prediction.x-tlo.x,2)/(2*pow(s_x, 2)) + (pow(associated_prediction.y-tlo.y,2)/(2*pow(s_y, 2))) ) );


            p.weight *= obs_w;
        }
    }


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;
    vector<Particle> new_particles;

    uniform_int_distribution<int> uniintdist(0, particles.size()-1);
    auto index = uniintdist(gen);


//    for(auto const & p : particles) {
//        cout << p.weight << endl;
//    }
     double max_weight = (*max_element(particles.begin(), particles.end(),
                                      [](Particle a, Particle b) {return a.weight <b.weight;})).weight;

     assert(max_weight > 0 );

     uniform_real_distribution<double> unirealdist(0.0, max_weight);

     double beta = 0.0;


     for (int i = 0; i < num_particles; i++) {
       beta += unirealdist(gen) * 2.0;
       while (beta > particles[index].weight) {
         beta -= particles[index].weight;
         index = (index + 1) % num_particles;
       }
       auto p = particles[index];
       p.id = i;
       new_particles.push_back(p);
     }

   particles = new_particles;

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
