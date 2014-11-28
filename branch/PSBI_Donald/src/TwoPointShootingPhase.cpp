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

    TwoPointShootingPhase::TwoPointShootingPhase(const int& j, const int& p, const missionoptions* options) :
        phase(j, p, options)
    {
        this->match_point_state.resize(7);
    }

    TwoPointShootingPhase::~TwoPointShootingPhase()
    {
        // default destructor is never used (I think it is superceded by the daughter class destructors)
    }

    void TwoPointShootingPhase::calcbounds_step_distribution_scale_factor(const string& prefix,
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
        if (options->step_size_distribution == 3 || options->step_size_distribution == 4)
        {
            Xlowerbounds->push_back(1.0);
            Xupperbounds->push_back(options->step_size_stdv_or_scale);
            Xdescriptions->push_back(prefix + "step size distribution standard deviation or scale factor");
        }
    }

    void TwoPointShootingPhase::calcbounds_LT_controls(const string& prefix,
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
        if (options->control_coordinate_system == 0) //Cartesian control
        {
            for (int w = 0; w < options->num_timesteps; ++w)
            {
                stringstream stepstream;
                stepstream << w;
                //u_x
                Xlowerbounds->push_back(-1.0);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_x");

                //u_y
                Xlowerbounds->push_back(-1.0);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_y");

                //u_z
                Xlowerbounds->push_back(-1.0);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_z");

                //for variable specific impulse propulsion systems, we must also encode the Isp for this time step
                if (options->engine_type == 4 || options->engine_type == 13)
                {
                    Xlowerbounds->push_back(options->IspLT_minimum);
                    Xupperbounds->push_back(options->IspLT);
                    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                }
                else if (options->engine_type == 12)
                {
                    Xlowerbounds->push_back(3000.0);
                    Xupperbounds->push_back(5000.0);
                    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                }

                //and the throttle magnitude constraint
                //throttle = 0
                Flowerbounds->push_back(0.0);
                Fupperbounds->push_back(1.0);
                Fdescriptions->push_back(prefix + "step " + stepstream.str() + " throttle magnitude constraint");

                //Jacobian entries for the throttle magnitude constraint
                //The throttle magnitude constraint is dependent only on the preceding three throttle components
                vector<int> step_G_indices;
                vector<double> dummyvalues(3);
                int vary_Isp_flag = (options->engine_type == 4 || options->engine_type == 12 || options->engine_type == 13) ? 1 : 0;
                for (size_t entry = Xdescriptions->size() - 1 - vary_Isp_flag; entry > Xdescriptions->size() - 4 - vary_Isp_flag; --entry)
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix + "step " << w << " throttle magnitude constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());

                    //store the position in the G vector
                    step_G_indices.push_back(iGfun->size() - 1);
                }
                control_vector_G_indices.push_back(step_G_indices);
            }
        }
        else if (options->control_coordinate_system == 1) //Polar control
        {
            for (int w = 0; w < options->num_timesteps; ++w)
            {
                stringstream stepstream;
                stepstream << w;
                //u_x
                Xlowerbounds->push_back(1.0e-10);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_Throttle");

                //u_y
                Xlowerbounds->push_back(0.0);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_u");

                //u_z
                Xlowerbounds->push_back(0.0);
                Xupperbounds->push_back(1.0);
                Xdescriptions->push_back(prefix + "step " + stepstream.str() + " u_v");

                //for variable specific impulse propulsion systems, we must also encode the Isp for this time step
                if (options->engine_type == 4 || options->engine_type == 13)
                {
                    Xlowerbounds->push_back(options->IspLT_minimum);
                    Xupperbounds->push_back(options->IspLT);
                    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                }
                else if (options->engine_type == 12)
                {
                    Xlowerbounds->push_back(3000.0);
                    Xupperbounds->push_back(5000.0);
                    Xdescriptions->push_back(prefix + "step " + stepstream.str() + " Isp");
                }
            }
        }
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
        
        if (p == 0 && options->allow_initial_mass_to_vary)
            this->G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier.resize(7);

        if (options->journey_variable_mass_increment[j])
            this->G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier.resize(7);

        if (options->objective_type == 13)
            this->G_index_of_derivative_of_match_point_with_respect_to_BOL_power.resize(7);

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
            //every match point constraint depends on all variables in the phase
            for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
            {
                //do NOT pick up entries with respect to flight time
                if (!((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024))
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());
                }

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
            //derivative with respect to times and epochs
            //the LAST entry is always the current phase flight time
            vector<int> state_G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables;
            vector<double> state_X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables;
            for (int entry = 0; entry <= Xdescriptions->size() - 1; ++entry)
            {
                if ((*Xdescriptions)[entry].find("time") < 1024 || (*Xdescriptions)[entry].find("epoch") < 1024)
                {
                    iGfun->push_back(Fdescriptions->size() - 1);
                    jGvar->push_back(entry);
                    stringstream EntryNameStream;
                    EntryNameStream << "Derivative of " << prefix << " patch point " << statename[state] << " constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                    Gdescriptions->push_back(EntryNameStream.str());
                    state_G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables.push_back(Gdescriptions->size() - 1);
                    state_X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                }
            }
            this->G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables.push_back(state_G_index_of_derivative_of_match_point_with_respect_to_flight_time_variables);
            this->X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables.push_back(state_X_scale_range_of_derivative_of_match_point_with_respect_to_flight_time_variables);

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
                                G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[state] = Gdescriptions->size() - 1;
                            }
                        }
                    }
                }
            }//end if (options->journey_variable_mass_increment[j])

            //check for derivatives with respect to escape and capture spirals in this journey and previous journeys
            this->find_dependencies_due_to_escape_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
            this->find_dependencies_due_to_capture_spiral(Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds, Xdescriptions, Fdescriptions, iAfun, jAvar, iGfun, jGvar, Adescriptions, Gdescriptions, j, p, options);
        }//end loop over state variables

        //Construct a helper array "rectangular prism" of G indices for match point derivatives with respect to state variables within the phase
        vector < vector <int> > timeslice;
        vector <int> scanline;
        vector<string> constraint_type;
        constraint_type.push_back(" patch point x constraint");
        constraint_type.push_back(" patch point y constraint");
        constraint_type.push_back(" patch point z constraint");
        constraint_type.push_back(" patch point xdot constraint");
        constraint_type.push_back(" patch point ydot constraint");
        constraint_type.push_back(" patch point zdot constraint");
        constraint_type.push_back(" patch point m constraint");

        for (int step = 0; step < options->num_timesteps + 2; ++step)
        {
            timeslice.clear();

            //loop over constraint type dimension
            for (int c = 0; c < 7; ++c)
            {
                scanline.clear();

                string constraintname = prefix + constraint_type[c];

                for (size_t entry = 0; entry < Gdescriptions->size(); ++entry)
                {
                    if ((*Gdescriptions)[entry].find(constraintname) < 1024)
                    {
                        if (step == 0) //derivatives with respect to initial velocity
                        {
                            //the first step of the first phase of a journey is abnormal because instead of encoding an XYZ vector, we have encoded a magnitude and two angles
                            if (p == 0)
                            {
                                if ((*Gdescriptions)[entry].find("magnitude of outgoing velocity asymptote") < 1024)
                                {
                                    scanline.push_back(entry);     //derivative with respect to magnitude
                                    scanline.push_back(entry + 1); //derivative with respect to RA
                                    scanline.push_back(entry + 2); //derivative with respect to DEC
                                }
                            }
                            //otherwise look for an XYZ initial velocity increment
                            else
                            {
                                if ((*Gdescriptions)[entry].find("initial velocity increment x") < 1024)
                                {
                                    scanline.push_back(entry);     //derivative with respect to initial velocity increment x
                                    scanline.push_back(entry + 1); //derivative with respect to initial velocity increment y
                                    scanline.push_back(entry + 2); //derivative with respect to initial velocity increment z
                                }
                            }


                            if ((*Gdescriptions)[entry].find("initial mass multiplier (0-1)") < 1024)
                            {
                                G_index_of_derivative_of_match_point_constraints_with_respect_to_mission_initial_mass_multiplier[c] = entry;
                            }
                            if ( /*p == 0 &&*/ (*Gdescriptions)[entry].find("journey initial mass scale factor") < 1024)
                            {
                                G_index_of_derivative_of_match_point_constraints_with_respect_to_journey_initial_mass_increment_multiplier[c] = entry;
                            }
                        }
                        //derivatives with respect to final velocity - these do not exist for a terminal rendezvous phase
                        else if (step == options->num_timesteps + 1 && (!(p == options->number_of_phases[j] - 1 && ((options->journey_arrival_type[j] == 1) || options->journey_arrival_type[j] == 3) || options->journey_arrival_type[j] == 5)))
                        {
                            if ((*Gdescriptions)[entry].find("terminal velocity increment x") < 1024)
                            {
                                scanline.push_back(entry);	   //derivative with respect to terminal velocity increment x
                                scanline.push_back(entry + 1); //derivative with respect to terminal velocity increment y
                                scanline.push_back(entry + 2); //derivative with respect to terminal velocity increment z
                            }
                        }
                        //derivatives with respect to the control parameters
                        else
                        {
                            stringstream stepstream;
                            stepstream << step - 1;

                            string controlname = prefix + "step " + stepstream.str() + " u_x";

                            if ((*Gdescriptions)[entry].find(controlname) < 1024)
                            {
                                scanline.push_back(entry);	   //derivative with respect to u_x
                                scanline.push_back(entry + 1); //derivative with respect to u_y
                                scanline.push_back(entry + 2); //derivative with respect to u_z
                            }
                        }


                    } //end if ( (*Gdescriptions)[entry].find(constraintname) )
                } //end loop over Gdescriptions

                timeslice.push_back(scanline);

            } //end loop over constraint entries
            match_point_constraint_G_indices.push_back(timeslice);

        } //end loop over time steps
    }

    void TwoPointShootingPhase::calcbounds_arrival_constraints(const string& prefix,
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
        if (p == options->number_of_phases[j] - 1)
        {
            if (options->journey_arrival_type[j] == 2) //intercept with bounded v-infinity
            {
                Flowerbounds->push_back((options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1])*(options->journey_final_velocity[j][0] / options->journey_final_velocity[j][1]) - 1);
                Fupperbounds->push_back(0.0);
                Fdescriptions->push_back(prefix + "arrival C3 constraint");

                //Jacobian entry for a bounded v-infinity intercept
                //this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        terminal_velocity_constraint_G_indices.push_back(iGfun->size() - 1);
                        terminal_velocity_constraint_X_indices.push_back(entry);
                        terminal_velocity_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
            }

            if (options->journey_arrival_declination_constraint_flag[j] && (options->journey_arrival_type[j] == 0 || options->journey_arrival_type[j] == 2)) //intercept with bounded v-infinity or orbit insertion
            {
                Flowerbounds->push_back(options->journey_arrival_declination_bounds[j][0] / options->journey_arrival_declination_bounds[j][1] - 1.0);
                Fupperbounds->push_back(0.0);
                Fdescriptions->push_back(prefix + "arrival declination constraint");

                //Jacobian entry for arrival declination constraint
                //this is a nonlinear constraint dependent only on the terminal velocity vector for this phase
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " arrival v-infinity constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        arrival_declination_constraint_G_indices.push_back(iGfun->size() - 1);
                        arrival_declination_constraint_X_indices.push_back(entry);
                        arrival_declination_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
            }

            else if (options->journey_arrival_type[j] == 6) //enforce escape constraint, E = 0
            {
                Flowerbounds->push_back(-math::SMALL);
                Fupperbounds->push_back(math::LARGE);
                Fdescriptions->push_back(prefix + "arrival escape condition (zero energy constraint)");

                //Jacobian entry for zero energy condition
                //this is a nonlinear constraint dependent on the terminal velocity vector
                //note that it is NOT dependent on position because these phases always end at the SOI and SOI distance is constant
                for (size_t entry = first_X_entry_in_phase; entry < Xdescriptions->size(); ++entry)
                {
                    if ((*Xdescriptions)[entry].find("terminal velocity") < 1024 || (*Xdescriptions)[entry].find("right boundary point") < 1024)
                    {
                        iGfun->push_back(Fdescriptions->size() - 1);
                        jGvar->push_back(entry);
                        stringstream EntryNameStream;
                        EntryNameStream << "Derivative of " << prefix << " escape constraint constraint F[" << Fdescriptions->size() - 1 << "] with respect to X[" << entry << "]: " << (*Xdescriptions)[entry];
                        Gdescriptions->push_back(EntryNameStream.str());
                        escape_constraint_G_indices.push_back(iGfun->size() - 1);
                        escape_constraint_X_indices.push_back(entry);
                        escape_constraint_X_scale_ranges.push_back((*Xupperbounds)[entry] - (*Xlowerbounds)[entry]);
                    }
                }
            }
        }
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

    void TwoPointShootingPhase::process_arrival(double* current_state,
                                                double* current_deltaV,
                                                double* boundary2_state,
                                                double* current_epoch,
                                                double* X,
                                                int* Xindex,
                                                double* F,
                                                int* Findex,
                                                double* G,
                                                const int& j,
                                                const int& p,
                                                const bool& needG,
                                                missionoptions* options,
                                                EMTG::Astrodynamics::universe* Universe)
    {
        if (options->journey_arrival_type[j] == 3 || options->journey_arrival_type[j] == 5)
            dV_arrival_magnitude = 0.0;

        //note that "3" for journey_arrival_type indicates a "low-thrust rendezvous," which means we are there already and we don't need to do anything
        else
        {
            //compute the arrival deltaV
            if (boundary2_location_code > 0) //ending at body
                dV_arrival_magnitude = phase::process_arrival(  state_at_end_of_phase + 3,
                                                                boundary2_state,
                                                                current_state + 3,
                                                                current_epoch,
                                                                Body2->mu,
                                                                Body2->r_SOI,
                                                                F,
                                                                Findex,
                                                                j,
                                                                options,
                                                                Universe);
            else //arriving at a boundary point in free space
                dV_arrival_magnitude = phase::process_arrival(  state_at_end_of_phase + 3,
                                                                boundary2_state,
                                                                current_state + 3,
                                                                current_epoch,
                                                                Universe->mu,
                                                                Universe->r_SOI,
                                                                F,
                                                                Findex,
                                                                j,
                                                                options,
                                                                Universe);

            //drop the electric propulsion stage
            state_at_end_of_phase[6] -= options->EP_dry_mass;

            //apply the arrival burn
            state_at_end_of_phase[6] *= exp(-dV_arrival_magnitude * 1000 / ((options->IspChem > 0 ? options->IspChem : 1e-6)* options->g0));
            *current_deltaV += dV_arrival_magnitude;

            //if this is a terminal intercept, we need to compute derivatives
            if (needG && options->journey_arrival_type[j] == 2)
            {
                double C3_desired = options->journey_final_velocity[j][1] * options->journey_final_velocity[j][1];
                //derivative with respect to x component of terminal velocity
                G[terminal_velocity_constraint_G_indices[0]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[0] * X[terminal_velocity_constraint_X_indices[0]] / C3_desired;

                //derivative with respect to y component of terminal velocity
                G[terminal_velocity_constraint_G_indices[1]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[1] * X[terminal_velocity_constraint_X_indices[1]] / C3_desired;

                //derivative with respect to z component of terminal velocity
                G[terminal_velocity_constraint_G_indices[2]] = 2.0 * terminal_velocity_constraint_X_scale_ranges[2] * X[terminal_velocity_constraint_X_indices[2]] / C3_desired;
            }
            else if (needG && options->journey_arrival_type[j] == 6)
            {
                double r = math::norm(boundary2_state, 3) / Universe->LU;
                double denominator = r * r * r * Universe->LU * Universe->LU;
                double LUTU2 = (Universe->TU / Universe->LU) * (Universe->TU / Universe->LU);
                double sqLU = sqrt(Universe->LU);

                //derivative with respect to x component of position
                G[escape_constraint_G_indices[0]] = escape_constraint_X_scale_ranges[0] * X[escape_constraint_X_indices[0]] / denominator;

                //derivative with respect to x component of position
                G[escape_constraint_G_indices[1]] = escape_constraint_X_scale_ranges[1] * X[escape_constraint_X_indices[1]] / denominator;

                //derivative with respect to x component of position
                G[escape_constraint_G_indices[2]] = escape_constraint_X_scale_ranges[2] * X[escape_constraint_X_indices[2]] / denominator;

                //derivative with respect to x component of terminal velocity
                G[escape_constraint_G_indices[3]] = escape_constraint_X_scale_ranges[3] * X[escape_constraint_X_indices[3]] * LUTU2;

                //derivative with respect to y component of terminal velocity
                G[escape_constraint_G_indices[4]] = escape_constraint_X_scale_ranges[4] * X[escape_constraint_X_indices[4]] * LUTU2;

                //derivative with respect to z component of terminal velocity
                G[escape_constraint_G_indices[5]] = escape_constraint_X_scale_ranges[5] * X[escape_constraint_X_indices[5]] * LUTU2;
            }
        }
    }
}//close namespace EMTG