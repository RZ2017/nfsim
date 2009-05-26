/*
 * CompartmentReaction.h
 *
 *  Created on: Apr 13, 2009
 *      Author: Kelly Stanton
 *
 *      The CompartmentReaction class specifies an extended version of the ReactionClass
 *      it adds functionality for storing molecules and the molecule lists in a list
 *      of compartments.  This is still an abstract class and will need to be derived
 *      for specific reaction types
 */


#ifndef COMPARTMENTREACTION_H_
#define COMPARTMENTREACTION_H_

#include "../NFreactions.hh"

#include "Compartment.h"
#include "CompartmentSubPropensity.h"
#include <vector>
using namespace NFcore;
using namespace std;

class CompartmentReactantList;
class CompartmentSubPropensity;
class CompartmentReaction: public NFcore::ReactionClass {
	friend class CompartmentReactantList;
	friend class CompartmentSubPropensity;
public:
	//Constructor
	CompartmentReaction(string name, double baseRate, TransformationSet *transformationSet);
	//Destructor
	virtual ~CompartmentReaction();
	//Initilization
	virtual void init();
	//Finalization
	virtual void prepareForSimulation();
	//Adds a molecule to this reaction rule
	virtual bool tryToAdd(Molecule *m, unsigned int reactantPos);
	//Removese a molecule from this reaction rule
	virtual void remove(Molecule *m, unsigned int reactantPos);
	//Gets the propensity for this reaction rule to fire
	virtual double update_a();

	virtual void notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex);
	//Gets the number of reactants
	virtual unsigned int getReactantCount(unsigned int reactantIndex) const;
	//prints information about this reaction.  Mainly used for debugging
	virtual void printFullDetails() const;

	//restrict this reaction to firing only in the specified compartment
	//this is done by changing the propensity of each compartment
	void restrictToCompartment(unsigned int compartmentId);

	// Changes the compartment the specified molecule is in for this reaction
	//void moveMolToCompartment(Molecule* m, unsigned int oldCompartmentId, unsigned int newCompartmentId, unsigned int reactantPos);

	//Must call before using the compartmentReaction class
	static void SetNumCompartments(unsigned int nCompartments);
	static void addConnectivity(unsigned int compartmentId1, unsigned int compartmentId2);


protected:

	//add an interaction between compartments
	void addCompartmentInteraction(unsigned int CompartmentIdList[]);

	//Choose a mappingSet at random
	virtual void pickMappingSets(double randNumber) const;

	//note maps are access in log(n) time
	map<unsigned int,CompartmentReactantList*> m_mapCompartmentList;

	//vector of interactions that this reaction can handle
	vector<CompartmentSubPropensity*> m_vectInteractionList;

	//Number of compartments in the system
	static unsigned int nCompartments;
	//keeps track of compartment connections
	static map<unsigned int,vector<unsigned int> > m_mapCompartmentConnectivity;
};

#endif /* COMPARTMENTREACTION_H_ */



