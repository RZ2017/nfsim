#include "reaction.hh"
#include "../../NFutil/setting.hh" //razi added for debugging purpose, last update 2017-3-29

#define DEBUG_MESSAGE 0//was 0  Razi changed to show more debugging messages


using namespace std;
using namespace NFcore;


#ifdef RHS_FUNC
//should also accept list of local functions and list of PointerNames for each of the functions...
RHSRxnClass::RHSRxnClass(   //Razi: this version supports RHS functions
		string name,
		double baseRate,
		string baseRateName,
		bool checkProducts,
		TransformationSet *transformationSet,
		CompositeFunction *function,
		vector <string> &lfArgumentPointerNameList, System *s) :
	ReactionClass(name,baseRate,baseRateName,checkProducts,transformationSet,s)
{
	bool verbose=false;
	if (DEBUG_ACTIVE & CREATE_REACTION)
		verbose = system->getverbose();

	if (verbose) {cout<<"\n\tRHS RXN:"<<name<<" with RHS composite functions:"<< function->getName() << "and "<< n_reactants <<" reactants is created.\n"; mypause(-1);}

	this->reactionType = RHS_RXN;
	this->checkProducts = checkProducts;
	reactantLists = new ReactantList *[n_reactants];
	check_mappingSet = new MappingSet *[n_reactants];

	//Set up the reactantLists
	for(unsigned int r=0; r<n_reactants; r++)
		reactantLists[r]=(new ReactantList(r,transformationSet,25));
	

	if (verbose) cout<<"RHSreaction: I developed everything up to here. More may be needed."; //mypause(-1); //exit(0);

	//Set the actual function
	this->cfo = function;
	cfo->RHS = true;
	//Initialize a to zero
	this->a=0; //later check

}




RHSRxnClass::~RHSRxnClass() {

	for(unsigned int r=0; r<n_reactants; r++) {
		delete reactantLists[r];
	}

	delete [] reactantLists;

/*
	delete [] argIndexIntoMappingSet;
	delete [] argMappedMolecule;
	delete [] argScope;
*/

}

void RHSRxnClass::init() {

	//Here we have to tell the molecules that they are part of this function
	//and for single molecule functions, we have to tell them also that they are in
	//this function, so they need to update thier value should they be transformed
	for(unsigned int r=0; r<n_reactants; r++)
	{
		reactantTemplates[r]->getMoleculeType()->addReactionClass(this,r);
	}
}




void RHSRxnClass::remove(Molecule *m, unsigned int reactantPos)
{

	cerr<<"irrelevant function remove for RHSRxnClass\n";
	//First a bit of error checking...
	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	{
		cout<<"Error removing molecule from a reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
		exit(1);
	}


	//Get the specified reactantList
	ReactantList *rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	bool isInRxn = (m->getRxnListMappingId(rxnIndex)>=0);


	if(isInRxn)
	{
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}
}





int RHSRxnClass::checkForEquality(Molecule *m, MappingSet* ms, int rxnIndex, ReactantList* reactantList){
	//Razi: Comes from BASIC RXN, check if any modifications are required
	/*
	Check if mapping set clashes with any of the mapping sets already in reactantList
	*/
	
	set<int> tempSet = m->getRxnListMappingSet(rxnIndex);
	for(set<int>::iterator it= tempSet.begin();it!= tempSet.end(); ++it){
		MappingSet* ms2 = reactantList->getMappingSet(*it);
		if(MappingSet::checkForEquality(ms,ms2)){
			return *it;
		}
	}
	return -1;
}




bool RHSRxnClass::tryToAdd(Molecule *m, unsigned int reactantPos) {
	//based on Basic RXN
	//First a bit of error checking, that you should skip unless we are debugging...
	//	if(reactantPos<0 || reactantPos>=n_reactants || m==NULL)
	//	{
	//		cout<<"Error adding molecule to reaction!!  Invalid molecule or reactant position given.  Quitting."<<endl;
	//		exit(1);
	//	}

	//Get the specified reactantList
	rl = reactantLists[reactantPos];

	//Check if the molecule is in this list
	int rxnIndex = m->getMoleculeType()->getRxnIndex(this,reactantPos);
	//cout<<" got mappingSetId: " << m->getRxnListMappingId(rxnIndex)<<" size: " <<rl->size()<<endl;
	//cout<< " testing whether to add molecule ";
	//m->printDetails();
	//cout<<" ... as a normal reaction "<<this->name<<endl;


	//If this reaction has multiple instances, we always remove them all!
	// then we remap because other mappings may have changed.  Yes, this may
	// be more ineffecient, but it is the fast implementation
	if(rl->getHasClonedMappings()) {
		while(m->getRxnListMappingId(rxnIndex)>=0) {
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
			//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}
	}

	/*
	//Here we get the standard update...
	if(m->getRxnListMappingId(rxnIndex)>=0) //If we are in this reaction...
	{
		if(!reactantTemplates[reactantPos]->compare(m)) {
			//	cout<<"Removing molecule "<<m->getUniqueID()<<" which was at mappingSet: "<<m->getRxnListMappingId(rxnIndex)<<endl;
			rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
			m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
		}

	} else {
		//Try to map it!
		ms = rl->pushNextAvailableMappingSet();
		if(!reactantTemplates[reactantPos]->compare(m,rl,ms)) {
			//we must remove, if we did not match.  This will also remove
			//everything that was cloned off of the mapping set
			rl->removeMappingSet(ms->getId());
		} else {
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}
	}
	*/

	//Here we get the standard update...
	while(m->getRxnListMappingId(rxnIndex)>=0) {
		rl->removeMappingSet(m->getRxnListMappingId(rxnIndex));
		m->deleteRxnListMappingId(rxnIndex,m->getRxnListMappingId(rxnIndex));
		//m->setRxnListMappingId(rxnIndex,Molecule::NOT_IN_RXN);
	}


	//Try to map it!
	ms = rl->pushNextAvailableMappingSet();
	symmetricMappingSet.clear();
	comparisonResult = reactantTemplates[reactantPos]->compare(m,rl,ms,false,&symmetricMappingSet);
	if(!comparisonResult) {
		//cout << "no mapping in normal reaction, remove"<<endl;
		//we must remove, if we did not match.  This will also remove
		//everything that was cloned off of the mapping set
		rl->removeMappingSet(ms->getId());
		//JJT: removes any symmetric mapping sets that might have been added since we are not using them
		for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
			rl->removeMappingSet((*it)->getId());
		}
	} else {
		//cout << "should be in normal reaction, confirm push"<<endl;
		//ms->printDetails();

		//TODO: it is necessary to remove elements that are not used anymore from the rl as well as from the m
		//for that
		//m->setRxnListMappingId(rxnIndex,-1);

		if (symmetricMappingSet.size() > 0){
            rl->removeMappingSet(ms->getId());
			for(vector<MappingSet *>::iterator it=symmetricMappingSet.begin();it!=symmetricMappingSet.end();++it){
					//XXX: JJT this is a band-aid, symmetricMappingSet should not have repeated elements in the first place
					int mapIndex = checkForEquality(m,*it,rxnIndex,rl);
					if(mapIndex >= 0){
						rl->removeMappingSet((*it)->getId());
					}
					else{
						m->setRxnListMappingId(rxnIndex,(*it)->getId());
					}
            }
		}
		else{
			m->setRxnListMappingId(rxnIndex,ms->getId());
		}

	}

	return true;
}





int RHSRxnClass::getReactantCount(unsigned int reactantIndex) const
{
	//Razi: Check if this is needed
	return isPopulationType[reactantIndex] ?
		       reactantLists[reactantIndex]->getPopulation()
	         : reactantLists[reactantIndex]->size();
}


int RHSRxnClass::getCorrectedReactantCount(unsigned int reactantIndex) const
{
	//Razi: Check if this is needed
	return isPopulationType[reactantIndex] ?
			   std::max( reactantLists[reactantIndex]->getPopulation()
			             - identicalPopCountCorrection[reactantIndex], 0 )
			 : reactantLists[reactantIndex]->size();
}


//This function takes a given mappingset and looks up the value of its local
//functions based on the local functions that were defined
double RHSRxnClass::evaluateLocalFunctions(MappingSet *ms)
{
	cerr<< "RHSRxnClass::evaluateLocalFunctions not developed yet.\n";
	return -1;
	//Razi: This may need substantial changes: Since we will need to apply the function on the reaction products
	// after firing a reaction, [and not on the input reactants]


	//Go through each function, and set the value of the function
	//this->argMappedMolecule
	//cout<<"\t\t\t\tRHSRxnClass::evaluateLocalFunctions()"<<endl;
	//cout<<"dor is reevaluating its function."<<endl;

	//Grab the molecules needed for the local function to evaluate
	for(int i=0; i<this->n_argMolecules; i++) {
		//cout<<"here."<<endl;
		//cout<<"\t\t\t\t\t"<<i<<": argMappedMolecule="<<argMappedMolecule[i]<<" argIndexIntoMappingSet="<<argIndexIntoMappingSet[i]<<endl;
		this->argMappedMolecule[i] = ms->get(this->argIndexIntoMappingSet[i])->getMolecule();
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeType="<<argMappedMolecule[i]->getMoleculeTypeName()<<endl;
		//cout<<"\t\t\t\t\t"<<"argMappedMoleculeScope="<<argScope[i]<<endl;
	}

	//cout<<"done setting molecules, so know calling the composite function evaluate method."<<endl;
	int * reactantCounts = new int[this->n_reactants];
	for(unsigned int r=0; r<n_reactants; r++) {
		reactantCounts[r]=reactantLists[r]->size();
	}

	double value = this->cfo->evaluateOn(argMappedMolecule,argScope, reactantCounts, n_reactants);
	delete [] reactantCounts;
	//cout<<"\t\t\t\t\t"<<"composite function value="<<value<<endl;

	return value;

	/*Molecule

	for(int i=0; i<(signed)lfList.size(); i++) {
		Molecule *molObject = ms->get(this->indexIntoMappingSet.at(i))->getMolecule();
		int index = lfList.at(i)->getIndexOfTypeIFunctionValue(molObject);
		this->localFunctionValue.at(i)=molObject->getLocalFunctionValue(index);
		//cout<<"found that local function: "<<getName()<<" evaluates to: " <<localFunctionValue.at(i)<<endl;
	}
	return this->localFunctionValue.at(0);
	*/
}


double RHSRxnClass::update_a() {
	bool verbose = false;
	if (res == 0)
		res = cfo->evaluateOnRHS(0, 0, CompositeFunction::EvalConditionalPart, false); 
if (DEBUG_ACTIVE & SHOW_SIM)
	verbose = system->getverbose();
#ifdef FIX_A   //Razi: Propensity of a reaction should be perhaps the minimum of reactant counts and not multiplication of it, later check
	double c;
	a = 1.0;
	a = getCorrectedReactantCount(0);
	if (n_reactants >= 1){
		for(unsigned int i=1; i<n_reactants; i++) {
			c=getCorrectedReactantCount(i);	a=(a<c)?a:c;    //min(a,getCorrectedReactantCount(i)
			cout<<"\t update a for RHS reaction:"<< this->getName()<<"  i:"<<i<<" count:"<<c<< "  a:min(a,c):"<<a<<endl;
		}
	}
	a*=baseRate;
	cout<<"\t updating finished: a for reaction after applying base rate is:"<< a<<endl;

#else

	if(this->totalRateFlag) {
		a=baseRate;
		for(unsigned int i=0; i<n_reactants; i++)
			if(getCorrectedReactantCount(i)==0) a=0.0;

	// Use the standard microscopic rate
	} else {
		a = 1.0;
		for(unsigned int i=0; i<n_reactants; i++) {
			a*=getCorrectedReactantCount(i);
		}
		a*=baseRate;
		if(res!=0)
			a*=res;
	}

#endif


	if (verbose){// && (this->getName().compare("XbPlusYbBind")==0) && 0){ //don't show anymore
		int i,j,k;
		cout <<"\tupdate_a() is called for RHS reaction:"<<this->name << " totalRateFlag:"<<totalRateFlag<<"   Base rate:"<<baseRate<<endl;
		for(i=0; i<n_reactants; i++){
			j = getCorrectedReactantCount(i);
			cout<<"\tReactantlist "<<i<<"  Population Type: "<< isPopulationType[i] <<" Corrected count: "<< j<<endl;
		}
		cout <<"\tFinal a is :"<<a<<endl<<endl;		mypause(500);
	}
	return a;
}

void RHSRxnClass::pickMappingSets(double randNumber) const
{
	//Razi: Source is BasicRxnClass::pickMappingSets check for modifications

	//Note here that we completely ignore the argument.  The argument is only
	//used for DOR reactions because we need that number to select the reactant to fire
	for(unsigned int i=0; i<n_reactants; i++)
	{
		if ( isPopulationType[i] ) {
			reactantLists[i]->pickRandomFromPopulation(mappingSet[i]);
		} else {
			reactantLists[i]->pickRandom(mappingSet[i]);
		}
	}
}

void RHSRxnClass::notifyRateFactorChange(Molecule * m, int reactantIndex, int rxnListIndex) {
	cerr<<"notifyRateFactorChange is not developed or is an invalid function for RHS reaction!!!\n."; exit(1);
}





void RHSRxnClass::printDetails() const
{
	cout<<"RHSRxnClass: " << name <<"  ( baseRate="<<baseRate<<",  a="<<a<<", fired="<<fireCounter<<" times )"<<endl;
	for(unsigned int r=0; r<n_reactants; r++)
	{
		cout<<"      -|"<< this->getReactantCount(r)<<" mappings|\t";
		cout<<this->reactantTemplates[r]->getPatternString()<<"\n";
	}

/*
	if (DEBUG_ACTIVE & SHOW_FIRE) //Razi: print details of ractants: for some reasons it does not work, later check (e.g.null pointers)
		for(unsigned int i=0; i<n_reactants; i++){
			try{
				reactantLists[i]->printDetails();
			}catch(...){
				cout<<"Error occurred while printing details of the DOR reaction.\n";
			}
	}
	*/

	if(n_reactants==0)
		cout<<"      >No Reactants: so this rule either creates new species or does nothing."<<endl;
}



#endif
