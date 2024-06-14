
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
#include "TMVA/CrossValidation.h"
using namespace TMVA;

void TMVARegression( TString RECO_files, Int_t species )
{
    // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
    // if you use your private .rootrc, or run from a different directory, please copy the
    // corresponding lines from .rootrc
    // methods to be processed can be given as an argument; use format:
    //
    //     mylinux~> root -l TMVARegression.CUndefined control sequence \"
    //
    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDEFoam"]         = 0;
    Use["KNN"]             = 0;
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    //
    // Neural Network
    Use["MLP"]             = 1;
    Use["DNN"]             = 0;
    //
    // Support Vector Machine
    Use["SVM"]             = 0;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 0;
    Use["BDTG"]            = 0;
    // ---------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "==> Start TMVARegression" << std::endl;
    // Select methods (don't look at this code - not of interest)
    // if (myMethodList != "") {
    //     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    //     std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    //     for (UInt_t i=0; i<mlist.size(); i++) {
    //         std::string regMethod(mlist[i]);
    //         if (Use.find(regMethod) == Use.end()) {
    //         std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
    //         for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
    //         std::cout << std::endl;
    //         return;
    //         }
    //         Use[regMethod] = 1;
    //     }
    // }
    // --------------------------------------------------------------------------------------------------

    Bool_t isKtautau = false;
    Bool_t isDpDmK = false;
    Bool_t isD0D0K = false;

    if( (species == 10) || (species == 11) || (species == 12) || (species == 1) || (species == 2) || (species == 3) )
    {
        isKtautau = true;
    }
    else if( (species == 4) || (species == 5) || (species == 6) )
    {
        isDpDmK = true;
    }
    else if( (species == 9) || (species == 0) || (species == -1) )
    {
        isD0D0K = true;
    }

    // Here the preparation phase begins
    // Create a new root output file
    TFile* outputFile = new TFile( Form("/panfs/felician/B2Ktautau/workflow/make_MLP_regression/TMVAReg_%i.root", species), "RECREATE");

    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory will
    // then run the performance analysis for you.
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weight/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in
    // front of the "Silent" argument in the option string

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
    // If you wish to modify default settings
    // (please check "src/Config.h" to see all available global options)
    //
    //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
    //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    if(isD0D0K)
    {
        dataloader->AddVariable("df_mprime_1", "PVx", 'F');
        dataloader->AddVariable("df_mprime_2", "PVy", 'F');
        dataloader->AddVariable("df_mprime_3", "PVz", 'F');
        dataloader->AddVariable("df_mprime_4", "DV1x", 'F');
        dataloader->AddVariable("df_mprime_5", "DV1y", 'F');
        dataloader->AddVariable("df_mprime_6", "DV1z", 'F');
        dataloader->AddVariable("df_mprime_7", "p3pi1x", 'F');
        dataloader->AddVariable("df_mprime_8", "p3pi1y", 'F');
        dataloader->AddVariable("df_mprime_9", "p3pi1z", 'F');
        dataloader->AddVariable("df_mprime_10", "E3pi1", 'F');
        dataloader->AddVariable("df_mprime_11", "DV2x", 'F');
        dataloader->AddVariable("df_mprime_12", "DV2y", 'F');
        dataloader->AddVariable("df_mprime_13", "DV2z", 'F');
        dataloader->AddVariable("df_mprime_14", "p3pi2x", 'F');
        dataloader->AddVariable("df_mprime_15", "p3pi2y", 'F');
        dataloader->AddVariable("df_mprime_16", "p3pi2z", 'F');
        dataloader->AddVariable("df_mprime_17", "E3pi2", 'F');
        dataloader->AddVariable("df_mprime_18", "RPx", 'F');
        dataloader->AddVariable("df_mprime_19", "RPy", 'F');
        dataloader->AddVariable("df_mprime_20", "pKx", 'F');
        dataloader->AddVariable("df_mprime_21", "pKy", 'F');
        dataloader->AddVariable("df_mprime_22", "pKz", 'F');
    }
    else
    {
        dataloader->AddVariable("df_m_1", "PVx", 'F');
        dataloader->AddVariable("df_m_2", "PVy", 'F');
        dataloader->AddVariable("df_m_3", "PVz", 'F');
        dataloader->AddVariable("df_m_4", "DV1x", 'F');
        dataloader->AddVariable("df_m_5", "DV1y", 'F');
        dataloader->AddVariable("df_m_6", "DV1z", 'F');
        dataloader->AddVariable("df_m_7", "p3pi1x", 'F');
        dataloader->AddVariable("df_m_8", "p3pi1y", 'F');
        dataloader->AddVariable("df_m_9", "p3pi1z", 'F');
        dataloader->AddVariable("df_m_10", "E3pi1", 'F');
        dataloader->AddVariable("df_m_11", "DV2x", 'F');
        dataloader->AddVariable("df_m_12", "DV2y", 'F');
        dataloader->AddVariable("df_m_13", "DV2z", 'F');
        dataloader->AddVariable("df_m_14", "p3pi2x", 'F');
        dataloader->AddVariable("df_m_15", "p3pi2y", 'F');
        dataloader->AddVariable("df_m_16", "p3pi2z", 'F');
        dataloader->AddVariable("df_m_17", "E3pi2", 'F');
        dataloader->AddVariable("df_m_18", "RPx", 'F');
        dataloader->AddVariable("df_m_19", "RPy", 'F');
        dataloader->AddVariable("df_m_20", "pKx", 'F');
        dataloader->AddVariable("df_m_21", "pKy", 'F');
        dataloader->AddVariable("df_m_22", "pKz", 'F');
    }

    dataloader->AddVariable("Kp_RP_Z", "RPz", 'F');
    dataloader->AddSpectator( "eventNumber",  "eventNumber", "", 'I' );

    // dataloader->AddVariable( "var1", "Variable 1", "units", 'F' );
    // dataloader->AddVariable( "var2", "Variable 2", "units", 'F' );
    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    // dataloader->AddSpectator( "spec1:=var1*2",  "Spectator 1", "units", 'F' );
    // dataloader->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );
    // Add the variable carrying the regression target
    if(isKtautau)
    {
        dataloader->AddTarget("taup_TRUEP_X");
        dataloader->AddTarget("taup_TRUEP_Y");
        dataloader->AddTarget("taup_TRUEP_Z");
        dataloader->AddTarget("taup_TRUEP_E");
        dataloader->AddTarget("taum_TRUEP_X");
        dataloader->AddTarget("taum_TRUEP_Y");
        dataloader->AddTarget("taum_TRUEP_Z");
        dataloader->AddTarget("taum_TRUEP_E");
    }
    else if(isDpDmK)
    {
        dataloader->AddTarget("Dp_TRUEP_X");
        dataloader->AddTarget("Dp_TRUEP_Y");
        dataloader->AddTarget("Dp_TRUEP_Z");
        dataloader->AddTarget("Dp_TRUEP_E");
        dataloader->AddTarget("Dm_TRUEP_X");
        dataloader->AddTarget("Dm_TRUEP_Y");
        dataloader->AddTarget("Dm_TRUEP_Z");
        dataloader->AddTarget("Dm_TRUEP_E");
    }
    else if(isD0D0K)
    {
        dataloader->AddTarget("D0_TRUEP_X");
        dataloader->AddTarget("D0_TRUEP_Y");
        dataloader->AddTarget("D0_TRUEP_Z");
        dataloader->AddTarget("D0_TRUEP_E");
        dataloader->AddTarget("D0bar_TRUEP_X");
        dataloader->AddTarget("D0bar_TRUEP_Y");
        dataloader->AddTarget("D0bar_TRUEP_Z");
        dataloader->AddTarget("D0bar_TRUEP_E");
    }

    // dataloader->AddTarget( targetVar );
    // It is also possible to declare additional targets for multi-dimensional regression, ie:
    //     factory->AddTarget( "fvalue2" );
    // BUT: this is currently ONLY implemented for MLP
    // Read training and test data (see TMVAClassification for reading ASCII files)
    // load the signal and background event samples from ROOT trees

    TFileCollection* fc = new TFileCollection("fc", "fc", RECO_files);
    TChain* regTree = new TChain("DecayTree");
    regTree->AddFileInfoList((TCollection*)fc->GetList());

    // global event weights per tree (see below for setting event-wise weights)
    Double_t regWeight  = 1.0;
    // You can add an arbitrary number of regression trees
    dataloader->AddRegressionTree( regTree, regWeight );
    // This would set individual event weights (the variables defined in the
    // expression need to exist in the original TTree)
    // dataloader->SetWeightExpression( "var1", "Regression" );
    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
    // tell the DataLoader to use all remaining events in the trees after training for testing:
    // dataloader->PrepareTrainingAndTestTree( mycut,
    //                                         "nTrain_Regression=0:"
    //                                         "nTest_Regression=0:"
    //                                         "NormMode=NumEvents:"
    //                                         "!V" );
    dataloader->PrepareTrainingAndTestTree( mycut,
                                        "nTest_Regression=1:"
                                        "NormMode=NumEvents:"
                                        "!V" );
    //
    //     dataloader->PrepareTrainingAndTestTree( mycut,
    //            "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
    // If no numbers of events are given, half of the events in the tree are used
    // for training, and the other half for testing:
    //
    //     dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
    // Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
    // PDE - RS method

    // TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
    //                                         "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );

    TString splitExpr = "int([eventNumber])\%int([numFolds])";

    TString opt = Form("!V:!Silent"
                       ":Color:DrawProgressBar"
                       ":AnalysisType=Regression"
                       ":NumFolds=2"
                       ":SplitType=Deterministic"
                       ":FoldFileOutput=False"
                       ":SplitExpr=%s", splitExpr.Data());
    TMVA::CrossValidation cv {Form("TMVARegression_%i", species), dataloader, outputFile, opt};
    // if (Use["PDERS"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kPDERS, "PDERS",
    //                         "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
    // // And the options strings for the MinMax and RMS methods, respectively:
    // //
    // //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    // //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    // if (Use["PDEFoam"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam",
    //             "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );
    // // K-Nearest Neighbour classifier (KNN)
    // if (Use["KNN"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kKNN, "KNN",
    //                         "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
    // // Linear discriminant
    // if (Use["LD"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kLD, "LD",
    //                         "!H:!V:VarTransform=None" );
    // // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    // if (Use["FDA_MC"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MC",
    //                         "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );
    // if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
    //     factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GA",
    //                         "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );
    // if (Use["FDA_MT"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MT",
    //                         "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
    // if (Use["FDA_GAMT"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GAMT",
    //                         "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
    // Neural network (MLP)
    if (Use["MLP"])
    {
        cv.BookMethod( TMVA::Types::kMLP, "MLP", 
                       "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:"
                       "HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:"
                       "SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
    }
    // if (Use["DNN"])
    // {
    // /*
    //     TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
    //     TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
    //     TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
    //     TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
    //     TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
    //     TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
    //     TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");
    // */
    //     TString layoutString ("Layout=TANH|100,LINEAR");
    //     TString training0 ("LearningRate=1e-5,Momentum=0.5,Repetitions=1,ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2");
    //     TString training1 ("LearningRate=1e-5,Momentum=0.9,Repetitions=1,ConvergenceSteps=170,BatchSize=30,TestRepetitions=7,WeightDecay=0.01,Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions=1");
    //     TString training2 ("LearningRate=1e-5,Momentum=0.3,Repetitions=1,ConvergenceSteps=150,BatchSize=40,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE");
    //     TString training3 ("LearningRate=1e-6,Momentum=0.1,Repetitions=1,ConvergenceSteps=500,BatchSize=100,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE");
    //     TString trainingStrategyString ("TrainingStrategy=");
    //     trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
    // //       TString trainingStrategyString ("TrainingStrategy=LearningRate=1e-1,Momentum=0.3,Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");
    //     TString nnOptions ("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
    // //       TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
    //     nnOptions.Append (":"); nnOptions.Append (layoutString);
    //     nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
    //     factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN", nnOptions ); // NN
    // }
    // // Support Vector Machine
    // if (Use["SVM"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
    // // Boosted Decision Trees
    // if (Use["BDT"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDT",
    //                         "!H:!V:NTrees=1000:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
    // if (Use["BDTG"])
    //     factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG",
    //                         "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
    // --------------------------------------------------------------------------------------------------
    // Now you can tell the factory to train, test, and evaluate the MVAs
    // Train MVAs using the set of training events
    // factory->TrainAllMethods();
    // // Evaluate all MVAs using the set of test events
    // factory->TestAllMethods();
    // // Evaluate and compare performance of all configured MVAs
    // factory->EvaluateAllMethods();

    cv.Evaluate();
    TMVA::CrossValidationResult results = cv.GetResults()[0];
    results.Print();
    // --------------------------------------------------------------
    // Save the output
    outputFile->Close();
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVARegression is done!" << std::endl;
    // delete factory;
    // delete dataloader;
    // Launch the GUI for the root macros
}