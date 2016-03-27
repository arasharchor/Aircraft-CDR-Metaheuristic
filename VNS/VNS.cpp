#define _USE_MATH_DEFINES

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <ctime>
#include "Aircraft.h"

using namespace std;

vector<Aircraft> VNS(
	vector<Aircraft> aircraft_vector,
	vector<vector<double> > neighborhood_vector,
	double expire_time,		// computing time limit;
	double improve_time,	// Computing time allowed since the last improvement is performed
	double ang,
	double penalty) {

	// Timer
	clock_t expire_start = clock();
	clock_t last_improvement;

	int iNieghborhoodCount = neighborhood_vector.size();
	int iAircraftCount = aircraft_vector.size();

	// Empty vector for angle variations
	vector<double> angle_vector(iAircraftCount, 0.0);

	// Objective function values
	double dOldObjective;
	double dNewObjective;

	vector<Aircraft> current_aircraft_vector = aircraft_vector;
	vector<Aircraft> new_aircraft_vector =
		FirstImprovement(aircraft_vector, ang, penalty);

	do {

		for (int ni = 0; ni < iNieghborhoodCount; ++ni) {

			vector<double> neighborhood = Shaking(neighborhood_vector.at(ni), ang);
			dOldObjective = ObjectiveFunction(current_aircraft_vector, angle_vector, penalty);
			new_aircraft_vector = FirstImprovement(aircraft_vector, ang, penalty);
			dNewObjective = ObjectiveFunction(new_aircraft_vector, neighborhood, penalty);

			// Make a move
			if (dNewObjective < dOldObjective) {
				current_aircraft_vector = new_aircraft_vector;
				ni = 0;
				last_improvement = clock();
			}
			else {
				continue;
			}

			if ((clock() - expire_start) / (double)CLOCKS_PER_SEC > improve_time)
				break;
		}

	} while (((clock() - expire_start) / (double)CLOCKS_PER_SEC) < expire_time);

	return current_aircraft_vector;
}

vector<double> Shaking(
	vector<double> neighborhood_vector,
	// Distance between two discretized values of ? (angle)
	double ang) {

	srand(time(NULL));
	double random = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

	int iNeighborCount = neighborhood_vector.size();
	int iIterationCount = iNeighborCount % 4;

	vector<double> temp_neighborhood_vector = neighborhood_vector;

	for (int ni = 0; ni < iIterationCount; ++ni) {
		if (random < 0.5)
			temp_neighborhood_vector.at(ni) += ang;
		else
			temp_neighborhood_vector.at(ni) -= ang;
	}

	return temp_neighborhood_vector;
}

vector<Aircraft> FirstImprovement(
	vector<Aircraft> aircraft_vector,
	// Distance between two discretized values of ? (angle)
	double ang,
	double penalty
	) {

	int iAircraftCount = aircraft_vector.size();
	vector<double> angle_vector(iAircraftCount, 0.0);
	vector<Aircraft> temp_aircraft_vector = aircraft_vector;

	// Objective function values
	double dOldObjective;
	double dNewObjective;

	for (int ai = 0; ai < iAircraftCount; ++ai) {

		dOldObjective = ObjectiveFunction(aircraft_vector, angle_vector, penalty);

		// Turn aircraft ai by ang degrees
		angle_vector.at(ai) = ang;

		dNewObjective = ObjectiveFunction(aircraft_vector, angle_vector, penalty);

		if (dNewObjective < dOldObjective) {
			// There's an improvement, continue with the next aircraft
			continue;
		}
		else {

			// Turn aircraft ai by -ang degrees
			angle_vector.at(ai) = -ang;
			dNewObjective = ObjectiveFunction(aircraft_vector, angle_vector, penalty);

			if (dNewObjective < dOldObjective) {
				// There's an improvement, continue with the next aircraft
				continue;
			}
			else {
				// Keep the current position
				angle_vector.at(ai) = 0.0;
			}
		}
	}

	for (int ai = 0; ai < iAircraftCount; ++ai) {
		temp_aircraft_vector.at(ai).dCurrentAngle += angle_vector.at(ai);
	}

	return temp_aircraft_vector;
}

double ObjectiveFunction(
	vector<Aircraft> aircraft_vector,	// Vector of aircrafts 
	vector<double> angle_vector,		// Vector of angle variations
	double penalty						// Penalty parameter of the infeasibility
	) {

	int iAircraftCount = aircraft_vector.size();
	int iAngleCount = angle_vector.size();

	// Aim is to minimize dObjectiveValue and maximize dPenaltyValue
	double dObjectiveValue = numeric_limits<double>::infinity();
	double dPenaltyValue = -numeric_limits<double>::infinity();

	double dLocalPenalty = 0.0;		// Penalty value for two aircrafts
	double dLocalObjective = 0.0;	// Objective value for two aircrafts

									// Ensure each aircraft has a corresponding angle

	assert(iAircraftCount == iAngleCount);

	for (int ai = 0; ai < iAircraftCount; ++ai) {

		dLocalPenalty = CalculateInfeasibilityForTwoAircrafts(
			aircraft_vector.at(ai),
			aircraft_vector.at((ai + 1) % iAircraftCount),
			angle_vector.at(ai),
			angle_vector.at((ai + 1) % iAngleCount));

		// Update penalty value
		if (dPenaltyValue < dLocalPenalty)
			dPenaltyValue = dLocalPenalty;

		dLocalObjective += CalculateFeasibilityForTwoAircrafts(
			aircraft_vector.at(ai),
			aircraft_vector.at((ai + 1) % iAircraftCount),
			angle_vector.at(ai),
			angle_vector.at((ai + 1) % iAngleCount));
	}

	return penalty * dPenaltyValue + dLocalObjective;
}

double CalculateFeasibilityForTwoAircrafts(
	Aircraft i,
	Aircraft j,
	double Mi,	// Angle variation for aircraft i
	double Mj	// Angle variation for aircraft j
	) {

	// There is a pathological situation, denomitaor equals to zero
	if (almost_equal((i.dVelocity * cos(i.dCurrentAngle + Mi) -
		j.dVelocity * cos(j.dCurrentAngle + Mj)), 0.0, 2)) {

		return	(i.dVelocity * sin(i.dCurrentAngle + Mi) -
			j.dVelocity * sin(j.dCurrentAngle + Mj)) /
			(i.dVelocity * cos(i.dCurrentAngle + Mi) -
				j.dVelocity * cos(j.dCurrentAngle + Mj));
	}
	else {

		return	(i.dVelocity * sin(i.dCurrentAngle + Mi + M_PI / 2) -
			j.dVelocity * sin(j.dCurrentAngle + Mj + M_PI / 2)) /
			(i.dVelocity * cos(i.dCurrentAngle + Mi + M_PI / 2) -
				j.dVelocity * cos(j.dCurrentAngle + Mj + M_PI / 2));
	}
}

double CalculateInfeasibilityForTwoAircrafts(
	Aircraft i,
	Aircraft j,
	double Mi,	// Angle variation for aircraft i
	double Mj	// Angle variation for aircraft j
	) {

	double dDistance = sqrt(pow(i.dPositionX - j.dPositionX, 2) +
		pow(i.dPositionY - j.dPositionY, 2));
	double dOmega = (i.dPositionY - j.dPositionY) / (i.dPositionX - j.dPositionX);
	double dAlpha = (i.dRadius + j.dRadius) / dDistance;
	double lij = dOmega + dAlpha;
	double gij = dOmega - dAlpha;
	double tij = 0.0;

	// There is a pathological situation, denomitaor equals to zero
	if (almost_equal((i.dVelocity * cos(i.dCurrentAngle + Mi) -
		j.dVelocity * cos(j.dCurrentAngle + Mj)), 0.0, 2)) {

		tij = (i.dVelocity * sin(i.dCurrentAngle + Mi + M_PI / 2) -
			j.dVelocity * sin(j.dCurrentAngle + Mj + M_PI / 2)) /
			(i.dVelocity * cos(i.dCurrentAngle + Mi + M_PI / 2) -
				j.dVelocity * cos(j.dCurrentAngle + Mj + M_PI / 2));
		return min(-1 / tan(lij) - tij, tij + 1 / tan(gij));
	}
	else {

		tij = (i.dVelocity * sin(i.dCurrentAngle + Mi) -
			j.dVelocity * sin(j.dCurrentAngle + Mj)) /
			(i.dVelocity * cos(i.dCurrentAngle + Mi) -
				j.dVelocity * cos(j.dCurrentAngle + Mj));
		return min(tan(lij) - tij, tij - tan(gij));
	}
}

// Reference: http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
	// the machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
		// unless the result is subnormal
		|| std::abs(x - y) < std::numeric_limits<T>::min();
}
