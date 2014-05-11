//EMTG outer-loop solution class
//for use with NSGAII and SGA
//M. Vavrina and J. Englander

#include "mission.h"
#include "EMTG_outerloop_solution.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "boost/algorithm/string.hpp"

namespace GeneticAlgorithm
{
	EMTG_outerloop_solution::EMTG_outerloop_solution(const std::vector<int>& X_in, const int& number_of_objectives)
	{
		this->X = X_in;
		this->fitness_values.resize(number_of_objectives, 1.0e+100);
	}

	bool EMTG_outerloop_solution::operator>(EMTG_outerloop_solution OtherSolution)
	{
		return this->compare_objective_greaterthan(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator<(EMTG_outerloop_solution OtherSolution)
	{
		return this->fitness_values[0] < OtherSolution.fitness_values[0];
	}

	bool EMTG_outerloop_solution::operator>=(EMTG_outerloop_solution OtherSolution)
	{
		return this->compare_objective_greaterthanorequalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator<=(EMTG_outerloop_solution OtherSolution)
	{
		return this->compare_objective_lessthanorequalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::operator==(EMTG_outerloop_solution OtherSolution)
	{
		return this->compare_objective_equalto(OtherSolution, 0, 1.0e-6);
	}

	bool EMTG_outerloop_solution::compare_objective_greaterthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] > tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_greaterthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] >= tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_lessthan(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] < tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_lessthanorequalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (this->fitness_values[index] - OtherSolution.fitness_values[index] <= tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_objective_equalto(EMTG_outerloop_solution& OtherSolution, size_t index, const double& tolfit)
	{
		return (fabs(this->fitness_values[index] - OtherSolution.fitness_values[index]) < tolfit ? true : false);
	}

	bool EMTG_outerloop_solution::compare_description(EMTG_outerloop_solution& OtherSolution)
	{
		return (this->description == OtherSolution.description ? true : false);
	}

	//function to evaluate a member of the population
	void EMTG_outerloop_solution::set_data_pointers(void* optionspointer, void* Universepointer)
	{
		options_base = (EMTG::missionoptions*) optionspointer;
		Universe = (boost::ptr_vector< EMTG::Astrodynamics::universe>*) Universepointer;
	}

	EMTG::missionoptions EMTG_outerloop_solution::parse_outer_loop_decision_vector()
	{
		//Step 1: initialize a counter and a description string
		int Xindex = 0;
		stringstream descriptionstream;
		descriptionstream << options_base->mission_name;

		//Step 2: create a new options structure as a copy of the original
		EMTG::missionoptions options = *options_base;
		options.run_outerloop = 0;
		options.quiet_basinhopping = options.quiet_outerloop;
		options.quiet_NLP = options.quiet_outerloop;

		//Step 3: decode decision vector and create an options structure
		//we must decode the decision vector in the order that we encoded it

		//global mission choices
		if (options.outerloop_vary_power && options.outerloop_power_choices.size() > 1)
		{
			options.power_at_1_AU = options.outerloop_power_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.power_at_1_AU << "kW";
		}
		if (options.outerloop_vary_launch_epoch && options.outerloop_launch_epoch_choices.size() > 1)
		{
			options.launch_window_open_date = options.outerloop_launch_epoch_choices[X[Xindex]] * 86400.0;
			++Xindex;
			descriptionstream << "_LD" << (int)(options.launch_window_open_date / 86400.0);
		}
		if (options.outerloop_vary_flight_time_upper_bound && options.outerloop_flight_time_upper_bound_choices.size() > 1)
		{
			options.total_flight_time_bounds[1] = options.outerloop_flight_time_upper_bound_choices[X[Xindex]] * 86400;

			if (X[Xindex] > 0 && options.outerloop_restrict_flight_time_lower_bound)
				options.total_flight_time_bounds[0] = options.outerloop_flight_time_upper_bound_choices[X[Xindex] - 1] * 86400.0;
			else
				options.total_flight_time_bounds[0] = 0.0;

			++Xindex;
			descriptionstream << "_" << (int)(options.total_flight_time_bounds[1] / 86400.0) << "d";
		}
		if (options.outerloop_vary_thruster_type && options.outerloop_thruster_type_choices.size() > 1)
		{
			options.engine_type = options.outerloop_thruster_type_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.thruster_names[options.engine_type];
		}
		if (options.outerloop_vary_number_of_thrusters && options.outerloop_number_of_thrusters_choices.size() > 1)
		{
			options.number_of_engines = options.outerloop_number_of_thrusters_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_nTh" << options.number_of_engines;
		}
		if (options.outerloop_vary_launch_vehicle && options.outerloop_launch_vehicle_choices.size() > 1)
		{
			options.LV_type = options.outerloop_launch_vehicle_choices[X[Xindex]];
			++Xindex;
			descriptionstream << "_" << options.LV_names[options.LV_type];
		}
		if (options.outerloop_vary_departure_C3 && options.outerloop_departure_C3_choices.size() > 1)
		{
			if (X[Xindex] > 0)
				options.journey_initial_impulse_bounds[0][0] = sqrt(options.outerloop_departure_C3_choices[X[Xindex] - 1]);
			else
				options.journey_initial_impulse_bounds[0][0] = 0.0;
			options.journey_initial_impulse_bounds[0][1] = sqrt(options.outerloop_departure_C3_choices[X[Xindex]]);

			descriptionstream << "_C3d" << options.outerloop_departure_C3_choices[X[Xindex]];
			++Xindex;
		}
		if (options.outerloop_vary_arrival_C3 && options.outerloop_arrival_C3_choices.size() > 1)
		{
			if (X[Xindex] > 0)
				options.journey_final_velocity.back()[0] = sqrt(options.outerloop_arrival_C3_choices[X[Xindex] - 1]);
			else
				options.journey_final_velocity.back()[0] = 0.0;
			options.journey_final_velocity.back()[1] = sqrt(options.outerloop_arrival_C3_choices[X[Xindex]]);

			descriptionstream << "_C3a" << options.outerloop_arrival_C3_choices[X[Xindex]];
			++Xindex;
		}
		//choices for each journey
		for (size_t j = 0; j < options.number_of_journeys; ++j)
		{
			vector<int> journey_sequence;
			descriptionstream << "_" << (*Universe)[j].central_body_name << "(";

			//journey destination
			if (options.outerloop_vary_journey_destination[j] && options.outerloop_journey_destination_choices[j].size() > 1)
			{
				//for successive journeys, starting location is the previous journey's destination
				if (j > 0)
				{
					options.destination_list[j][0] = options.destination_list[j - 1][1];
				}

				options.destination_list[j][1] = options.outerloop_journey_destination_choices[j][X[Xindex]];
				++Xindex;
			}

			//encode journey destination 1
			journey_sequence.push_back(options.destination_list[j][0]);
			switch (options.destination_list[j][0])
			{
			case -3: //begin at point on orbit
			{
						 descriptionstream << "o";
						 break;
			}
			case -2: //begin at fixed point
			{
						 descriptionstream << "p";
						 break;
			}
			case -1: //begin at SOI
			{
						 descriptionstream << "s";
						 break;
			}
			default:
				descriptionstream << (*Universe)[j].bodies[options.destination_list[j][0] - 1].short_name;
			}

			//flyby sequence using Englander and Conway "null gene" method
			if (options.outerloop_vary_journey_flyby_sequence[j] && options.outerloop_journey_flyby_sequence_choices[j].size() > 1)
			{
				for (size_t p = 0; p < options.outerloop_journey_maximum_number_of_flybys[j]; ++p)
				{
					//is this a valid flyby or has the outer-loop chosen a null code?
					if (X[Xindex] < options.outerloop_journey_flyby_sequence_choices[j].size())
					{
						journey_sequence.push_back(options.outerloop_journey_flyby_sequence_choices[j][X[Xindex]]);
						descriptionstream << (*Universe)[j].bodies[journey_sequence.back() - 1].short_name;
					}

					++Xindex;

				}
			}

			//encode journey destination 2
			journey_sequence.push_back(options.destination_list[j][1]);
			switch (options.destination_list[j][1])
			{
			case -3: //begin at point on orbit
			{
						 descriptionstream << "o";
						 break;
			}
			case -2: //begin at fixed point
			{
						 descriptionstream << "p";
						 break;
			}
			case -1: //begin at SOI
			{
						 descriptionstream << "s";
						 break;
			}
			default:
				descriptionstream << (*Universe)[j].bodies[options.destination_list[j][1] - 1].short_name;
			}

			descriptionstream << ")";

			options.sequence.push_back(journey_sequence);
			options.number_of_phases[j] = journey_sequence.size() - 1;
			vector<int> journey_phase_type(journey_sequence.size() - 1, options.mission_type);
			options.phase_type.push_back(journey_phase_type);
		}

		vector< vector<int> > sequence_input;
		vector< vector<int> > phase_type_input;

		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			vector<int> journey_sequence_input;
			vector<int> journey_phase_type_input;
			for (int p = 1; p < options.max_phases_per_journey + 1; ++p)
			{
				if (p < options.number_of_phases[j])
				{
					journey_sequence_input.push_back(options.sequence[j][p]);
					journey_phase_type_input.push_back(options.phase_type[j][p]);
				}
				else
				{
					journey_sequence_input.push_back(0);
					journey_phase_type_input.push_back(options.mission_type);
				}
			}
			sequence_input.push_back(journey_sequence_input);
			phase_type_input.push_back(journey_phase_type_input);
		}
		options.sequence_input.clear();
		options.phase_type_input.clear();
		options.sequence_input.push_back(sequence_input);
		options.number_of_trial_sequences = 1;
		options.phase_type_input = phase_type_input;

		options.mission_name = descriptionstream.str();
		options.description = "";
		this->description = options.mission_name;

		//if this is a previously evaluated solution that we want to take another cut at, seed the optimizer from its decision vector
		if (this->Xinner.size() > 0)
		{
			options.seed_MBH = 1;
			options.current_trialX = this->Xinner;
		}
		else
			options.seed_MBH = 0;

		return options;
	}

	void EMTG_outerloop_solution::evaluate()
	{
		//Step 1: create the options structure
		EMTG::missionoptions options = this->parse_outer_loop_decision_vector();

		//Step 2: instantiate a mission
		EMTG::mission TrialMission(&options, *Universe);

		//Step 3: run the inner-loop
		options.print_options_file(options.working_directory + "//" + this->description + ".emtgopt");
		if (!options.quiet_outerloop)
		{
			TrialMission.output_mission_tree(options.working_directory + "//" + this->description + ".emtgtree");
			cout << "Optimizing mission: " << TrialMission.options.description << endl;
		}

		//Step 3.1: Do a quick check - if any journey destination is visited twice (later make this a switch) then assign poor fitness values and return
		for (int j = 0; j < options.number_of_journeys; ++j)
		{
			for (int jj = 0; jj < options.number_of_journeys; ++jj)
			{
				if (!(jj == j) && (options.destination_list[j][1] == options.destination_list[jj][1]))
				{
					std::cout << "Mission " << TrialMission.options.description << " contains duplicate journeys. Discarding." << std::endl;
					this->innerloop_fitness = 1.0e+100;
					for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
						this->fitness_values[objective] = 1.0e+100;
					return;
				}
			}
		}
		//Step 3.2 if we made it this far then optimize
		TrialMission.optimize();

		//Step 4: extract the fitness values
		//if the trial mission has no feasible solution then set the fitness values very high
		if (TrialMission.number_of_solutions == 0 || TrialMission.options.outputfile.find("FAILURE") < 1024)
		{
			this->innerloop_fitness = 1.0e+100;
			for (size_t objective = 0; objective < options.outerloop_objective_function_choices.size(); ++objective)
				this->fitness_values[objective] = 1.0e+100;
			return;
		}
		else
		{
			this->Xinner = TrialMission.Xopt;
			for (int entry = 0; entry < TrialMission.Xdescriptions.size(); ++entry)
			{
				if (TrialMission.Xdescriptions[entry].find("epoch") < 1024 || TrialMission.Xdescriptions[entry].find("time") < 1024)
				{
					this->Xinner[entry] /= 86400.0;
				}
			}
			TrialMission.extract_objective_function_values(this->fitness_values);
			this->innerloop_fitness = TrialMission.F[0];
		}

		return;
	}

	void EMTG_outerloop_solution::write_csv_line(std::string filename)
	{
		std::ofstream outputfile(filename.c_str(), ios::out | ios::app);

		for (size_t gene = 0; gene < this->X.size(); ++gene)
			outputfile << this->X[gene] << ",";
		outputfile << this->description << "," << this->generation_found << "," << this->timestamp;
		for (size_t objective = 0; objective < this->fitness_values.size(); ++objective)
			outputfile << "," << this->fitness_values[objective];
		for (size_t Xinner_entry = 0; Xinner_entry < this->Xinner.size(); ++Xinner_entry)
			outputfile << "," << Xinner[Xinner_entry];
		outputfile << std::endl;

		outputfile.close();
	}

	void EMTG_outerloop_solution::clear_innerloop_solution()
	{
		this->Xinner.clear();
		this->description.clear();
		
		for (int obj = 0; obj < this->fitness_values.size(); ++obj)
			this->fitness_values[obj] = 1.0e+100;

		this->innerloop_fitness = 1.0e+100;
	}
}