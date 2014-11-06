//source file for TwoPointShootingPhase
//this is an intermediate phase type that one should never instantiate - it contains features common to MGALT, FBLT, MGANDSM (for now)

#include <iostream>
#include <fstream>
#include <sstream>

#include "TwoPointShootingPhase.h"
#include "Kepler_Lagrange_Laguerre_Conway_Der.h"
#include "EMTG_math.h"

namespace EMTG {
    TwoPointShootingPhase::TwoPointShootingPhase()
    {
        // default constructor is never used

    }

    TwoPointShootingPhase::~TwoPointShootingPhase()
    {
        // default destructor is never used (I think it is superceded by the daughter class destructors)
    }

    void TwoPointShootingPhase::calcbounds_LT_match_points(const string& prefix,
                                                            int first_X_entry_in_phase, 
                                                            vector<double>* Xupperbounds, 
                                                            vector<double>* Xlowerbounds, 
                                                            vector<double>* Fupperbounds, 
                                                            vector<double>* Flowerbounds, 
                                                            vector<string>* Xdescriptions,
                                                            vector<string>* Fdescriptions, 
                                                            vector<int>* iAfun, 
                                                            vector<int>* jAvar, 
                                                            vector<int>* iGfun, 
                                                            vector<int>* jGvar, 
                                                            vector<string>* Adescriptions,
                                                            vector<string>* Gdescriptions,
                                                            int j,
                                                            int p, 
                                                            EMTG::Astrodynamics::universe* Universe,
                                                            missionoptions* options)
    {
        G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass.resize(7);
        G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass.resize(7);

        vector<string> statename;
        statename.push_back("x");
        statename.push_back("y");
        statename.push_back("z");
        statename.push_back("xdot");
        statename.push_back("ydot");
        statename.push_back("zdot");
        statename.push_back("m");


        for (int state = 0; state < 7; ++state)
        {
            Flowerbounds->push_back(-math::SMALL);
            Fupperbounds->push_back(math::SMALL);
            Fdescriptions->push_back(prefix + "match point " + statename[state]);
            for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
            {
                iGfun->push_back(Fdescriptions->size() - 1);
                jGvar->push_back(entry);
                stringstream EntryNameStream;
                EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                Gdescriptions->push_back(EntryNameStream.str());

                if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
                {
                    G_index_of_derivative_of_match_point_constraints_with_respect_to_arrival_mass[state] = Gdescriptions->size() - 1;
                }
            }
            for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
            {
                if ((*Xdescriptions)[entry].find("arrival mass") < 1024)
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());

                    G_index_of_derivative_of_match_point_constraints_with_respect_to_initial_mass[state] = Gdescriptions->size() - 1;
                    break;
                }
            }
            //derivative with respect to times and epochs in previous journeys/phases
            for (int pj = 0; pj <= j; ++pj)
            {
                for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
                {
                    stringstream pprefix_stream;
                    pprefix_stream << "j" << pj << "p" << pp;
                    string pprefix = pprefix_stream.str();

                    for (int entry = 0; entry <= first_X_entry_in_phase; ++entry)
                    {
                        if ((*Xdescriptions)[entry].find(pprefix) < 1024 && ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024))
                        {
                            iGfun->push_back(Fdescriptions->size() - 1);
                            jGvar->push_back(entry);
                            stringstream EntryNameStream;
                            EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                            Gdescriptions->push_back(EntryNameStream.str());
                        }
                    }
                }
            }
            //all match point constraints have a dependency on the BOL power if it is a variable
            if (options->objective_type == 13)
            {
                if (state == 0)
                    G_index_of_derivative_of_match_point_with_respect_to_BOL_power.resize(7);

                for (size_t Xentry = 0; Xentry < Xdescriptions->size(); ++Xentry)
                {
                    if ((*Xdescriptions)[Xentry].find("engine input power (kW)") < 1024)
                    {
                        bool duplicateflag = false;
                        stringstream entry_tag_stream;
                        entry_tag_stream << "constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]";
                        for (int Gentry = Gdescriptions->size() - 1; Gentry >= 0; --Gentry)
                        {
                            if ((*Gdescriptions)[Gentry].find(entry_tag_stream.str()) < 1024)
                            {
                                duplicateflag = true;
                                G_index_of_derivative_of_match_point_with_respect_to_BOL_power[state] = Gentry;
                                this->power_range = (*Xupperbounds)[Xentry] - math::SMALL;
                                break;
                            }
                        }
                        if (!duplicateflag)
                        {
                            iGfun->push_back(Fdescriptions->size() - 1);
                            jGvar->push_back(Xentry);
                            stringstream EntryNameStream;
                            EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << Xentry << "]: " << (*Xdescriptions)[Xentry];
                            Gdescriptions->push_back(EntryNameStream.str());
                            G_index_of_derivative_of_match_point_with_respect_to_BOL_power[state] = Gdescriptions->size() - 1;
                            this->power_range = (*Xupperbounds)[Xentry] - math::SMALL;
                            break;
                        }
                    }
                }
            }
            //the mass constraint has a derivative with respect to all previous journey initial mass increment scale factors
            if (options->journey_variable_mass_increment[j])
            {
                for (int pj = 0; pj <= j; ++pj)
                {
                    for (int pp = 0; pp < (pj == j ? p : options->number_of_phases[pj]); ++pp)
                    {
                        stringstream pprefix_stream;
                        pprefix_stream << "j" << pj << "p" << pp;
                        string pprefix = pprefix_stream.str();

                        for (int entry = first_X_entry_in_phase; entry >= 0; --entry)
                        {
                            if ((*Xdescriptions)[entry].find(pprefix) < 1024 && (*Xdescriptions)[entry].find("journey initial mass scale factor") < 1024)
                            {
                                iGfun->push_back(Fdescriptions->size() - 1);
                                jGvar->push_back(entry);
                                stringstream EntryNameStream;
                                EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                                Gdescriptions->push_back(EntryNameStream.str());
                            }
                        }
                    }
                }
            }//end if (options->journey_variable_mass_increment[j])

            //check for derivatives with respect to escape and capture spirals in this journey and previous journeys
            this->find_dependencies_due_to_escape_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
            this->find_dependencies_due_to_capture_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
        }//end loop over state variables
    }

    void TwoPointShootingPhase::create_initial_guess(const int& desired_mission_type,
        const bool& VSI,
        double& current_epoch,
        const int& j,
        const int& p,
        vector<double>& NewX,
        int& NewXIndex,
        const vector<string>& NewXDescriptions,
        const missionoptions& options,
        const Astrodynamics::universe& Universe)
    {
        if (p == 0)
        {
            //the first entry is always the journey initial mass scale factor if applicable
            if (options.journey_starting_mass_increment[j] > 0.0 && options.journey_variable_mass_increment[j])
            {
                NewX.push_back(this->journey_initial_mass_increment_scale_factor);
                ++NewXIndex;
            }
            //unless this journey starts with a flyby
            if (!(options.journey_departure_type[j] == 3 || options.journey_departure_type[j] == 6))
            {
                //the next entry is always the launch date or stay time unless this journey starts with a flyby
                NewX.push_back((this->phase_start_epoch - current_epoch) / 86400.0);
                ++NewXIndex;

                //variable left boundary conditions
                if (options.destination_list[j][0] == -1)
                {
                    if (options.journey_departure_elements_type[j] == 0) //Cartesian
                    {
                        for (size_t state = 0; state < 6; ++state)
                        {
                            if (options.journey_departure_elements_vary_flag[j][state])
                            {
                                NewX.push_back(this->left_boundary_local_frame_state[state]);
                                ++NewXIndex;
                            }
                        }
                    }
                    else //COE
                    {
                        for (size_t state = 0; state < 6; ++state)
                        {
                            if (options.journey_departure_elements_vary_flag[j][state])
                            {
                                NewX.push_back(this->left_boundary_orbit_elements[state]);
                                ++NewXIndex;
                            }
                        }
                    }
                }
                //then encode the departure asymptote (three variables)
                if (!(options.journey_departure_type[j] == 5 || options.journey_departure_type[j] == 2))
                {
                    NewX.push_back(sqrt(this->C3_departure));
                    NewX.push_back(this->RA_departure);
                    NewX.push_back(this->DEC_departure);
                    NewXIndex += 3;

                    //mission initial mass multiplier if applicable
                    if (j == 0 && options.allow_initial_mass_to_vary)
                    {
                        NewX.push_back(this->mission_initial_mass_multiplier);
                        ++NewXIndex;
                    }
                }
            }
        }
        else
        {
            //encode flyby information
            if (desired_mission_type == 1)
            {
                //MGA-DSM phases want a periapse distance and b-plane angle - we're making up an angle for now but some day it should be calculated properly
                NewX.push_back(1.0 + this->flyby_altitude / this->Body1->radius);
                NewX.push_back(math::PI);
                NewXIndex += 2;
            }
            else if (desired_mission_type == 2 || desired_mission_type == 3 || desired_mission_type == 4 || desired_mission_type == 5)
            {
                //MGALT, FBLT, MGANDSM phases want an initial v-infinity vector
                NewX.push_back(this->V_infinity_out(0));
                NewX.push_back(this->V_infinity_out(1));
                NewX.push_back(this->V_infinity_out(2));
                NewXIndex += 3;
            }
        }

        //all phase types place the time of flight next
        NewX.push_back(this->TOF / 86400.0);
        ++NewXIndex;

        //variable right-hand boundary
        if (options.destination_list[j][1] == -1)
        {
            if (options.journey_arrival_elements_type[j] == 0) //Cartesian
            {
                for (size_t state = 0; state < 6; ++state)
                {
                    if (options.journey_arrival_elements_vary_flag[j][state])
                    {
                        NewX.push_back(this->right_boundary_local_frame_state[state]);
                        ++NewXIndex;
                    }
                }
            }
            else //COE
            {
                for (size_t state = 0; state < 6; ++state)
                {
                    if (options.journey_arrival_elements_vary_flag[j][state])
                    {
                        NewX.push_back(this->right_boundary_orbit_elements[state]);
                        ++NewXIndex;
                    }
                }
            }
        }

        //MGALT, FBLT, MGANDSM, and PSBI phases want the terminal velocity vector, if applicable, and then spacecraft mass at end of phase next
        if (desired_mission_type == 2 || desired_mission_type == 3 || desired_mission_type == 4 || desired_mission_type == 5)
        {
            if (p < options.number_of_phases[j] - 1 || options.journey_arrival_type[j] == 2)
            {
                double boundary_state[6];
                this->Body2->locate_body(this->phase_end_epoch, boundary_state, false, (missionoptions*)&options);
                NewX.push_back(this->V_infinity_out(0));
                NewX.push_back(this->V_infinity_out(1));
                NewX.push_back(this->V_infinity_out(2));
                NewXIndex += 3;
            }

            NewX.push_back(this->state_at_end_of_phase[6]);
            ++NewXIndex;
        }

        //MGADSM and MGANDSM phases want the burn index next
        //we'll pretend it's in the middle
        if (desired_mission_type == 1 || desired_mission_type == 4)
        {
            NewX.push_back(0.5);
            ++NewXIndex;
        }

        //MGALT, FBLT, and PSBI phases want the control values (and sometimes state values) next
        if (desired_mission_type == 2 || desired_mission_type == 3 || desired_mission_type == 5)
        {
            //set all of the control parameters to 0 for now
            //maybe later introduce some delta-v smoothing
            for (int step = 0; step < options.num_timesteps; ++step)
            {
                if (desired_mission_type == 5)
                {
                    double state_left[6];
                    if (options.mission_type == 2) //for MGALT we have to subtract off the maneuver
                    {
                        double unperturbed_state[6];
                        for (size_t state = 0; state < 3; ++state)
                        {
                            unperturbed_state[state] = this->spacecraft_state[step][state];
                            unperturbed_state[state + 3] = this->spacecraft_state[step][state + 3] - this->dV[step][state];
                        }

                        Kepler::Kepler_Lagrange_Laguerre_Conway_Der(unperturbed_state,
                            state_left,
                            Universe.mu,
                            Universe.LU,
                            -this->time_step_sizes[step] / 2.0);
                    }
                    for (size_t state = 0; state < 6; ++state)
                        NewX.push_back(state_left[state]);
                    NewX.push_back(this->spacecraft_state[step][6]);
                    NewXIndex += 7;
                }
                NewX.push_back(this->control[step][0]);
                NewX.push_back(this->control[step][1]);
                NewX.push_back(this->control[step][2]);
                NewXIndex += 3;
                //set Isp, if VSI, to options.IspLT
                if (VSI)
                {
                    NewX.push_back(this->available_Isp[step]);
                    ++NewXIndex;
                }
            }
        }

        //finally set the current epoch to the phase end epoch
        current_epoch = this->phase_end_epoch;
        return;
    }
}//close namespace EMTG