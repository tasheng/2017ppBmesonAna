//HltInfo
int           Bf_HLT_Run;
ULong64_t     Bf_HLT_Event;
int           Bf_HLT_LumiBlock;
void setHltBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Bf_HLT_Run);
  hltroot->SetBranchAddress("Event",&Bf_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Bf_HLT_LumiBlock);
}

//hiEvtInfo
unsigned int       Bf_HiTree_Run;
unsigned long long Bf_HiTree_Evt;
unsigned int       Bf_HiTree_Lumi;
int hiBin;

void setHiTreeBranch(TTree* hitreeroot)
{
  hitreeroot->SetBranchAddress("run",&Bf_HiTree_Run);
  hitreeroot->SetBranchAddress("evt",&Bf_HiTree_Evt);
  hitreeroot->SetBranchAddress("lumi",&Bf_HiTree_Lumi);
  hitreeroot->SetBranchAddress("hiBin",&hiBin);

}




bool pclusterCompatibilityFilter;
bool pprimaryVertexFilter;
bool phfCoincFilter2Th4;


void setSkimTreeBranch(TTree* skimTree)
{

	skimTree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
	skimTree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	skimTree->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);

}
