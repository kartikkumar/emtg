//source file for TwoPointShootingPhase
//this is an intermediate phase type that one should never instantiate - it contains features common to MGALT, FBLT, MGANDSM (for now)

#include <iostream>
#include <fstream>
#include <sstream>

#include "TwoPointShootingPhase.h"
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
}//close namespace EMTG