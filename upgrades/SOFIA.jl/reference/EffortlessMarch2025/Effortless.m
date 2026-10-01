(* ::Package:: *)

BeginPackage["Effortless`"]


(*Unprotect@@Names["Effortless`*"];
ClearAll@@Names["Effortless`*"];*)


RunEffortless::usage = "RunEffortless[SqrtList,EvenAlpha,SqrtRepl,vars] - constructs all possible odd letters for a given list of square roots SqrtList and a given list of even letters EvenAlpha. SqrtRepl contains replacement rules sending the abstract square roots in SqrtList to their expressions in terms of the variables in vars."
RunEffortlessIndependent::usage="RunEffortlessIndependent[SqrtList,EvenAlpha,SqrtRepl,vars] - constructs all possible independent odd letters for a given list of square roots SqrtList and a given list of even letters EvenAlpha. SqrtRepl contains replacement rules sending the abstract square roots in SqrtList to their expressions in terms of the variables in vars."
RunEffortlessIndependentFF::usage="RunEffortlessIndependentFF[SqrtList,EvenAlpha,SqrtRepl,vars] - uses FiniteFlow to construct all possible independent odd letters for a given list of square roots SqrtList and a given list of even letters EvenAlpha. SqrtRepl contains replacement rules sending the abstract square roots in SqrtList to their expressions in terms of the variables in vars."


FindOddLetters::usage = "FindOddLetters[CurrentSqrt,AllowedEvenAlpha,SqrtRepl,vars] - finds all odd letters constructed using the squareroot Sqrt which factorize in terms of the even alphabet EvenAlpha."
AllowedEvenLetters::usage = "AllowedEvenLetters[CurrentSqrt,EvenAlpha,SqrtRepl,vars] - finds all even letters in EvenAlpha which can possibly show up in a factorization of an odd letter involving Sqrt."
(*SortAlphabet::usage = "sorts the letters in EvenAlpha with respect to their polynomial order in the variables vars and returns the sorted alphabet as well as a map to the original alphabet."*)
(*SortAlphabetSublist::usage = "SortAlphabetSublist[EvenAlpha,vars] - sorts the letters in EvenAlpha with respect to their polynomial order in the variables vars and returns the sorted alphabet split in sublists corresponding to the polynomial order."*)
(*AlphabetPolynomialOrderKey::usage = "AlphabetPolynomialOrderKey[EvenAlpha,vars] - creates a key assigning every alphabet letter to its polynomial order in the variables vars."*)
ConstructProducts::usage = "ConstructProducts[CurrentSqrt,AllowedEvenAlpha,SqrtRepl,vars] - constructs all possible products of the right polynomial order for the square root Sqrt from the even subalphabet AllowedEvenAlpha."
AllowedProductQ::usage = "AllowedProductQ[CurrentSqrt,Prod,EvenAlpha,SqrtRepl,vars] - checks whether adding the abstract letter product Prod turns CurrentSqrt into a perfect square once EvenAlpha has been substituted."
AllowedLetterQ::usage = "AllowedLetterQ[SqrtArgument,EvenLetter,vars] - checks whether the letter EvenLetter is a letter that can possibly show up in the factorization of an odd letter involving SqrtArgument."
FindAllZeroLocations::usage = "FindAllZeroLocations[EvenAlpha] - finds the zero locations of all letters in EvenAlpha. (Only for letters which are polynomials of degree less than 5.)"
FactorOdd::usage="FactorOdd[OddLetter, SqrtRepl] - factorizes the product of numerator and denominator of OddLetter. SqrtRepl contains replacement rules for the appearing square root."
TakeDerivativeOfLetter::usage="TakeDerivativeOfLetter[Letter,SqrtRepl,vars,var] - takes the derivative of the letter Letter with respect to var. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
TakeDerivativeOfLetterModSqrt::usage="TakeDerivativeOfLetterModSqrt[Letter,SqrtRepl,vars,var] - takes the derivative of the letter Letter with respect to var. If Letter is an odd letter with sqrt Q, the derivative is multiplied by Q. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
FindIndependentLetters::usage="FindIndependentLetters[NewLetters,SqrtRepl,vars] - constructs the list of independent letters among NewLetters. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
FindIndependentLettersFF::usage="FindIndependentLettersFF[NewLetters,SqrtRepl,vars] - uses FiniteFlow to construct the list of independent letters among NewLetters. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
DecomposeInAlphabet::usage="DecomposeInAlphabet[AlphabetRepl,LetterToDecompose,SqrtRepl,vars] - decomposes the letter LetterToDecompose in terms of the alphabet provided in AlphabetRepl. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
DecomposeInAlphabetFF::usage="DecomposeInAlphabetFF[AlphabetRepl,LetterToDecompose,SqrtRepl,vars] - uses FiniteFlow to decompose the letter LetterToDecompose in terms of the alphabet provided in AlphabetRepl. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
FactorOddInAlphabet::usage="FactorOddInAlphabet[AlphabetRepl,OddLetter,SqrtRepl,vars] - factorizes the product of denominator and numerator of an odd letter OddLetter in terms of the alphabet provided in AlphabetRepl. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
FactorOddInAlphabetFF::usage="FactorOddInAlphabetFF[AlphabetRepl,OddLetter,SqrtRepl,vars] - uses FiniteFlow to factorize the product of denominator and numerator of an odd letter OddLetter in terms of the alphabet provided in AlphabetRepl. SqrtRepl contains replacement rules for the appearing square roots and vars is the list of variables."
(*TakeDerivativeOfOddModSqrt::usage=""*)
TakeDerivativeOfOdd::usage="TakeDerivativeOfOdd[Letter,SqrtRepl,vars,x] - takes the derivative of an odd letter Letter with respect to variable x."


ELImpossible::usage="Error Message for impossible decompositions."


Begin["`Private`"]


Options[RunEffortless]:={"ConstantCoefficientList"->{4},"CheckNumberProducts"->4,"CheckNumberLetters"->4,"Range"->1000,"Verbose"->True,"DoubleSquareRoots"->False,"ExtraPolynomials"->{1}}
RunEffortless[SqrtList_,EvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
	Module[{AllowedLetters,OddLetters,AllOddLetters={},StartTime=AbsoluteTime[],LocalTime,AllZeroLocations,sqrtList},
			AllZeroLocations=FindAllZeroLocations[EvenAlpha];
			If[OptionValue["DoubleSquareRoots"]==True,sqrtList=Join[SqrtList,Table[SqrtList[[i]]SqrtList[[j]],{i,1,Length[SqrtList]},{j,i+1,Length[SqrtList]}]//Flatten],sqrtList=SqrtList];
			Do[
				Print["For coefficient c = "<>ToString[Currentc]<>":"];
				Do[
					LocalTime=AbsoluteTime[];
					If[OptionValue["Verbose"]==True,Print["--------- Analyzing square root: "<>ToString[CurrentSqrt]]];
	
					OddLetters=Flatten[Table[
										AllowedLetters=DeleteCases[AllowedEvenLetters[CurrentPoly CurrentSqrt,EvenAlpha,SqrtRepl,vars,"CheckNumberLetters"->OptionValue["CheckNumberLetters"],"AllZeroLocations"->AllZeroLocations,"Range"->OptionValue["Range"]],Rule[_,Log[CurrentPoly]]|Rule[_,Log[-CurrentPoly]],All];
										If[OptionValue["Verbose"]==True,Print["Found "<>ToString[Length[Flatten[AllowedLetters]]]<>" allowed letters for odd letters involving "<>ToString[CurrentPoly CurrentSqrt]<>"."]];
										FindOddLetters[CurrentPoly CurrentSqrt,AllowedLetters,SqrtRepl,vars,"ConstantCoefficient"->Currentc,"CheckNumberProducts"->OptionValue["CheckNumberProducts"],"Range"->OptionValue["Range"]],
									{CurrentPoly,OptionValue["ExtraPolynomials"]}]];
					If[OptionValue["Verbose"]==True,Print["Constructed "<>ToString[Length[OddLetters]]<>" odd letters in "<>ToString[AbsoluteTime[]-LocalTime]<>" seconds."]];
					AppendTo[AllOddLetters,OddLetters];
					,{CurrentSqrt,sqrtList}
				];
				,{Currentc,OptionValue["ConstantCoefficientList"]}
			];
			
				Print["Constructed "<>ToString[Length[Flatten[AllOddLetters]]]<>" odd letters, total runtime: "<>ToString[AbsoluteTime[]-StartTime]<>" seconds."];
				AllOddLetters
			];



CheckInputEvenAlphaQ[EvenAlpha_]:=
Module[{},
If[Head[EvenAlpha]===List,If[Head[EvenAlpha[[1]]]===Rule,{True},{False,"ERROR: Even Alphabet is not a list of rules."}]
,{False,"ERROR: Even Alphabet is not a list."}]	
]


Options[RunEffortlessIndependent]:={"ConstantCoefficientList"->{4},"CheckNumberProducts"->4,"CheckNumberLetters"->4,"Range"->1000,"NumberOfNumericConfig"->50,"Verbose"->True,"DoubleSquareRoots"->False,"ExtraPolynomials"->{1}}
RunEffortlessIndependent[SqrtList_,EvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
Module[{oddList,independentoddList,CEvenAlpha=CheckInputEvenAlphaQ[EvenAlpha],IndependentOddLetters,sqrtList},
	If[!CEvenAlpha[[1]],
		CEvenAlpha[[2]],
		oddList=RunEffortless[SqrtList,EvenAlpha,SqrtRepl,vars,"ConstantCoefficientList"->OptionValue["ConstantCoefficientList"],"CheckNumberProducts"->OptionValue["CheckNumberProducts"],"CheckNumberLetters"->OptionValue["CheckNumberLetters"],"Range"->OptionValue["Range"],"Verbose"->OptionValue["Verbose"],"DoubleSquareRoots"->OptionValue["DoubleSquareRoots"],"ExtraPolynomials"->OptionValue["ExtraPolynomials"]];
		If[OptionValue["DoubleSquareRoots"]==True,sqrtList=Join[SqrtList,Table[SqrtList[[i]]SqrtList[[j]],{i,1,Length[SqrtList]},{j,i+1,Length[SqrtList]}]//Flatten],sqrtList=SqrtList];
		oddList=Table[oddList[[j]]//Table[Head[EvenAlpha[[1,1]]][sqrtList[[j]]][i]->#[[i]],{i,1,Length[#]}]&,{j,1,Length[oddList]}]//DeleteCases[#,{}]&;
		If[Length[Flatten[oddList]]==0,IndependentOddLetters={};,
		independentoddList=Table[FindIndependentLetters[oddList[[i]],SqrtRepl,vars,"NumberOfNumericConfig"->OptionValue["NumberOfNumericConfig"],"Verbose"->OptionValue["Verbose"],"CurrentSqrt"->(Head[oddList[[i,1,1]]]/.{Head[EvenAlpha[[1,1]]][smth_]->smth})],{i,1,Length[oddList]}];
		IndependentOddLetters=Flatten[independentoddList]/.Rule->List//Transpose//Last;];
		Print["Constructed "<>ToString[Length[Flatten[IndependentOddLetters]]]<>" independent odd letters."];
		IndependentOddLetters
		]
	]


Options[RunEffortlessIndependentFF]:={"ConstantCoefficientList"->{4},"CheckNumberProducts"->4,"CheckNumberLetters"->4,"Range"->1000,"Verbose"->True,"DoubleSquareRoots"->False,"ExtraPolynomials"->{1}}
RunEffortlessIndependentFF[SqrtList_,EvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=Module[{oddList,independentoddList,IndependentOddLetters,sqrtList},
	oddList=RunEffortless[SqrtList,EvenAlpha,SqrtRepl,vars,"ConstantCoefficientList"->OptionValue["ConstantCoefficientList"],"CheckNumberProducts"->OptionValue["CheckNumberProducts"],"CheckNumberLetters"->OptionValue["CheckNumberLetters"],"Range"->OptionValue["Range"],"Verbose"->OptionValue["Verbose"],"DoubleSquareRoots"->OptionValue["DoubleSquareRoots"],"ExtraPolynomials"->OptionValue["ExtraPolynomials"]];
	If[OptionValue["DoubleSquareRoots"]==True,sqrtList=Join[SqrtList,Table[SqrtList[[i]]SqrtList[[j]],{i,1,Length[SqrtList]},{j,i+1,Length[SqrtList]}]//Flatten],sqrtList=SqrtList];
	oddList=Table[oddList[[j]]//Table[Head[EvenAlpha[[1,1]]][sqrtList[[j]]][i]->#[[i]],{i,1,Length[#]}]&,{j,1,Length[oddList]}]//DeleteCases[#,{}]&;
	If[Length[Flatten[oddList]]==0,IndependentOddLetters={};,
		independentoddList=Table[FindIndependentLettersFF[oddList[[i]],SqrtRepl,vars,"Verbose"->OptionValue["Verbose"],"CurrentSqrt"->(Head[oddList[[i,1,1]]]/.{Head[EvenAlpha[[1,1]]][smth_]->smth})],{i,1,Length[oddList]}];
		IndependentOddLetters=Flatten[independentoddList]/.Rule->List//Transpose//Last;];
	(*independentoddList=Table[FindIndependentLettersFF[oddList[[i]],SqrtRepl,vars,"Verbose"->OptionValue["Verbose"],"CurrentSqrt"->(Head[oddList[[i,1,1]]]/.{Head[EvenAlpha[[1,1]]][smth_]->smth})],{i,1,Length[oddList]}];
	IndependentOddLetters=Flatten[independentoddList]/.Rule->List//Transpose//Last;*)
	Print["Constructed "<>ToString[Length[Flatten[IndependentOddLetters]]]<>" independent odd letters."];
	IndependentOddLetters
		]


Options[FindOddLetters]:={"ConstantCoefficient" -> 4,"CheckNumberProducts" -> 4,"Range"->1000};
FindOddLetters[CurrentSqrt_,AllowedEvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
	Module[{ProdCoeff,odd,Poly},
		ProdCoeff=ConstructProducts[CurrentSqrt,AllowedEvenAlpha,SqrtRepl,vars,"ConstantCoefficient"->OptionValue["ConstantCoefficient"]];
		odd={};
		Do[
			If[AllowedProductQ[CurrentSqrt,CurrentProd,AllowedEvenAlpha,SqrtRepl,vars,"CheckNumberProducts"->OptionValue["CheckNumberProducts"],"Range"->OptionValue["Range"]],
				Poly=Sqrt[Factor[(CurrentSqrt)^2+CurrentProd/.SqrtRepl/.Flatten[AllowedEvenAlpha]]]//PowerExpand;
				AppendTo[odd,Log[(Poly-CurrentSqrt)/(Poly+CurrentSqrt)]];
				]
		,{CurrentProd,ProdCoeff}];
	odd];


Options[AllowedProductQ]:={"CheckNumberProducts" -> 4,"Range"->1000};
AllowedProductQ[CurrentSqrt_,Prod_,EvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
	Module[{i,IntegerPoints,QQ=(CurrentSqrt)^2/.SqrtRepl,ProdExpl=Prod/.Flatten[EvenAlpha],point},
		IntegerPoints=Table[vars[[i]]->RandomInteger[{-OptionValue["Range"],OptionValue["Range"]}],{j,1,OptionValue["CheckNumberProducts"]},{i,1,Length[vars]}];
		And@@Table[IntegerQ[Sqrt[QQ+ProdExpl]/.point],{point,IntegerPoints}]
		];


Options[ConstructProducts]:={"ConstantCoefficient" -> 4};
ConstructProducts[CurrentSqrt_,AllowedEvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
	If[Flatten[AllowedEvenAlpha]==={},{},
		Module[{\[Alpha]=Flatten[AllowedEvenAlpha][[1,1,0]],qPartitions,q,prodlist,Prod,PList,red,j},
			q=Max[Cases[CoefficientRules[(CurrentSqrt)^2/.SqrtRepl],v_?VectorQ:>Total[v],2]];
			Table[red[j]=AllowedEvenAlpha[[j]][[;;,1,1]],{j,1,Length[AllowedEvenAlpha]}];
			qPartitions=IntegerPartitions[q,All,Table[j,{j,1,Length[AllowedEvenAlpha]}]];
			prodlist=Table[Tuples[Table[red[entry],{entry,part}]],{part,qPartitions}]//Flatten[#,1]&;
			Prod=Times@@(Exp[\[Alpha][#]]&/@#)&/@prodlist;
			Join[OptionValue["ConstantCoefficient"] Prod,-OptionValue["ConstantCoefficient"] Prod]//DeleteDuplicates
			]
		];


FindAllZeroLocations[EvenAlpha_]:=
	Module[{},
			Table[entry[[1]]->If[((Max[Cases[CoefficientRules[Exp[entry[[2]]]],v_?VectorQ:>Total[v],2]]))<5,Solve[Exp[entry[[2]]]==0],{}],{entry,EvenAlpha}]
			];


Options[AllowedEvenLetters]:={"AllZeroLocations"->{},"CheckNumberLetters"->4,"Range"->1000};
AllowedEvenLetters[CurrentSqrt_,EvenAlpha_,SqrtRepl_,vars_,OptionsPattern[]]:=
	Module[{AllowedLetters={},QQ=CurrentSqrt/.SqrtRepl},
		Do[
			If[OptionValue["AllZeroLocations"]==={},
				If[AllowedLetterQ[QQ,Exp[entry[[2]]],vars,"CheckNumberLetters"->OptionValue["CheckNumberLetters"],"Range"->OptionValue["Range"]],
					AppendTo[AllowedLetters,entry];
					];,
				If[AllowedLetterQ[QQ,Exp[entry[[2]]],vars,"CheckNumberLetters"->OptionValue["CheckNumberLetters"],"Range"->OptionValue["Range"],"ZeroLocations"->entry[[1]]/.OptionValue["AllZeroLocations"]],
					AppendTo[AllowedLetters,entry];
					];];,	
			{entry,EvenAlpha}];
		If[AllowedLetters==={},{},SortAlphabetSublist[AllowedLetters,vars]]
		];


Options[AllowedLetterQ]:={"CheckNumberLetters"->4,"Range"->1000,"ZeroLocations"->{}}
AllowedLetterQ[QQ_,EvenLetter_,vars_,OptionsPattern[]]:=
	Module[{IntegerPoints,ZeroLocations,ValidTestQ},
		If[((Max[Cases[CoefficientRules[EvenLetter],v_?VectorQ:>Total[v],2]]))<5,
			IntegerPoints=Table[vars[[i]]->(-1)^RandomInteger[{1,2}]RandomInteger[{1,OptionValue["Range"]}],{j,1,OptionValue["CheckNumberLetters"]},{i,1,Length[vars]}];
			If[OptionValue["ZeroLocations"]==={},
				ZeroLocations=Solve[EvenLetter==0];,
				ZeroLocations=OptionValue["ZeroLocations"];
				];
			ValidTestQ=And@@#&/@(And[IntegerQ[Numerator[#]],IntegerQ[Denominator[#]]]&/@#&/@Flatten/@Transpose[Quiet[vars/.ZeroLocations/.IntegerPoints]]);
			ZeroLocations=Last/@DeleteCases[Transpose[{ValidTestQ,ZeroLocations}],{False,_}];
			And@@(Flatten[QQ/.ZeroLocations/.IntegerPoints]//And[IntegerQ[Numerator[#]],IntegerQ[Denominator[#]]]&/@#&),
			True]
		];


SortAlphabetSublist[EvenAlpha_,vars_]:=
	Module[{APOK=AlphabetPolynomialOrderKey[EvenAlpha,vars],MaxOrder,jj},
		MaxOrder=Max[APOK/.Rule->List//Transpose//Last];
		Table[Select[EvenAlpha,(#[[1]]/.APOK)==jj&],{jj,1,MaxOrder}]
		]


AlphabetPolynomialOrderKey[EvenAlpha_,vars_]:=
	Module[{jj},
		EvenAlpha//Table[#[[jj]]//(#[[1]]->(Max[Cases[CoefficientRules[Exp[#[[2]]]],v_?VectorQ:>Total[v],2]]))&,{jj,1,Length[EvenAlpha]}]&
		]


Options[FindIndependentLetters]:={"NumberOfNumericConfig"->50,"Verbose"->True,"CurrentSqrt"->{},"Range"->100}
FindIndependentLetters[NewLetters_,SqrtRepl_,vars_,OptionsPattern[]]:=Module[{list1,list2,decomp},
		list1=NewLetters; (* list of input log letters (to be added) *)
		list2={NewLetters[[1]]}; (* output, list of independent letters (can also be a list of already independent letters from previous analysis) *)

		For[kk=2,kk<=Length[list1],kk++,
			decomp=DecomposeInAlphabet[list2,list1[[kk]][[2]],SqrtRepl,vars,"NumberOfNumericConfig"->OptionValue["NumberOfNumericConfig"],"CurrentSqrt"->OptionValue["CurrentSqrt"],"Range"->OptionValue["Range"]];
			If[decomp===ELImpossible,AppendTo[list2,list1[[kk]]];If[OptionValue["Verbose"]==True,Print["New letter found at position ",kk," : ",list1[[kk]]]]]
			];
	list2
	]


Options[FindIndependentLettersFF]:={"Verbose"->True,"CurrentSqrt"->{}}
FindIndependentLettersFF[NewLetters_,SqrtRepl_,vars_,OptionsPattern[]]:=Module[{list1,list2,decomp},
		list1=NewLetters; (* list of input log letters (to be added) *)
		list2={NewLetters[[1]]}; (* output, list of independent letters (can also be a list of already independent letters from previous analysis) *)

		For[kk=2,kk<=Length[list1],kk++,
			decomp=DecomposeInAlphabetFF[list2,list1[[kk]][[2]],SqrtRepl,vars,"CurrentSqrt"->OptionValue["CurrentSqrt"]];
			If[decomp===ELImpossible,AppendTo[list2,list1[[kk]]];If[OptionValue["Verbose"]==True,Print["New letter found at position ",kk," : ",list1[[kk]]]]]
			];
	list2
	]


Options[DecomposeInAlphabet]:={"NumberOfNumericConfig"->50,"CurrentSqrt"->{},"Range"->100}
DecomposeInAlphabet[AlphabetRepl_,LetterToDecompose_,SqrtRepl_,vars_,OptionsPattern[]]:=Module[{c,ansatzabstract,ansatz,Dansatz, DansatzNlist,DansatzN,NumericList,sol,solution,Alphabet=Last/@AlphabetRepl,\[Alpha]=Head[First[First[AlphabetRepl]]]},
	ansatzabstract=Sum[c[i]AlphabetRepl[[i,1]], {i,1,Length[AlphabetRepl]}];
	Dansatz[var_]:=Dansatz[var]=Sum[c[i]TakeDerivativeOfLetterModSqrt[Alphabet[[i]],SqrtRepl,vars,var,"CurrentSqrt"->OptionValue["CurrentSqrt"]],{i,1,Length[Alphabet]}];
	Table[Dansatz[var];, {var,vars}];
	Do[NumericList=Table[var->RandomReal[OptionValue["Range"]],OptionValue["NumberOfNumericConfig"],{var,vars}];
		DansatzNlist[var_]:=DansatzNlist[var]=Dansatz[var]/.NumericList;Table[DansatzNlist[var];, {var,vars}];
		DansatzN=Table[Dansatz[var], {var,vars}]/.NumericList;sol[i]=Table[TakeDerivativeOfLetterModSqrt[LetterToDecompose,SqrtRepl,vars,var,"CurrentSqrt"->OptionValue["CurrentSqrt"]], {var,vars}]/.NumericList//NSolve[(#)==(DansatzN)]&;
	,{i,1,10}];
	solution=Table[sol[i],{i,1,10}]//Chop//Rationalize//DeleteDuplicates;
	If[solution=={{}}, ELImpossible,ansatzabstract/.Flatten[solution]]&//First
		]


Options[DecomposeInAlphabetFF]:={"CurrentSqrt"->{}}
DecomposeInAlphabetFF[AlphabetRepl_,LetterToDecompose_,SqrtRepl_,vars_,OptionsPattern[]]:=Module[{del,c,cSol,ansatzabstract,ansatz,Dansatz,DansatzNlist,DansatzN,NumericList,Alphabet=Last/@AlphabetRepl,\[Alpha]=Head[First[First[AlphabetRepl]]],EqSys},
	If[Needs["FiniteFlow`"]==$Failed,Print["FiniteFlow not found."]];
	ansatzabstract=Sum[c[i]AlphabetRepl[[i,1]], {i,1,Length[AlphabetRepl]}];
	Dansatz[var_]:=Dansatz[var]=Sum[c[i]TakeDerivativeOfLetterModSqrt[Alphabet[[i]],SqrtRepl,vars,var,"CurrentSqrt"->OptionValue["CurrentSqrt"]],{i,1,Length[Alphabet]}];
	EqSys=Sum[del[var]TakeDerivativeOfLetterModSqrt[LetterToDecompose,SqrtRepl,vars,var,"CurrentSqrt"->OptionValue["CurrentSqrt"]], {var,vars}]//(#)-(Sum[del[var]Dansatz[var],{var,vars}])&;
	cSol=FiniteFlow`FFLinearFit[{},EqSys,0,Join[vars,Table[del[var],{var,vars}],{}],Table[c[i],{i,1,Length[AlphabetRepl]}]];
	If[cSol===FiniteFlow`FFImpossible, Return[ELImpossible],ansatzabstract/.cSol]
		]


FactorOdd[OddLetter_,SqrtRepl_]:=
	Module[{expLett=Exp[OddLetter]},
	Numerator[expLett]Denominator[expLett]//.SqrtRepl//Expand//Factor
	]


Options[FactorOddInAlphabet]:={"NumberOfNumericConfig"->50,"Range"->100}
FactorOddInAlphabet[AlphabetRepl_,OddLetter_,SqrtRepl_,vars_,OptionsPattern[]]:=
	Module[{factoredOdd},
	factoredOdd=FactorOdd[OddLetter,SqrtRepl];
	DecomposeInAlphabet[AlphabetRepl,Log[factoredOdd],SqrtRepl,vars,"NumberOfNumericConfig"->OptionValue["NumberOfNumericConfig"],"Range"->OptionValue["Range"]]//Variables
	]


FactorOddInAlphabetFF[AlphabetRepl_,OddLetter_,SqrtRepl_,vars_]:=
	Module[{factoredOdd},
	factoredOdd=FactorOdd[OddLetter,SqrtRepl];
	DecomposeInAlphabetFF[AlphabetRepl,Log[factoredOdd],SqrtRepl,vars]//Variables
	]


TakeDerivativeOfLetter[Letter_,SqrtRepl_,vars_,var_]:=
	Module[{expLett,Q},
		expLett=Letter//Cases[#,Log[__],All]&//First//Exp;
		Q=expLett//Numerator//#/.ToRules[vars==0]&;
		If[Q===0,
			D[Letter,var],
			TakeDerivativeOfOdd[Letter,SqrtRepl,vars,var]
			]
		]


Options[TakeDerivativeOfLetterModSqrt]:={"CurrentSqrt"->{}}
TakeDerivativeOfLetterModSqrt[Letter_,SqrtRepl_,vars_,var_,OptionsPattern[]]:=
	Module[{expLett,Q},
		expLett=Letter//Exp;
		If[OptionValue["CurrentSqrt"]==={},Q=expLett//Numerator//#/.ToRules[vars==0]/.{-smth_->smth}&,Q=OptionValue["CurrentSqrt"]];
		If[Q===0,
			D[Letter,var],
			TakeDerivativeOfOddModSqrt[Letter,SqrtRepl,vars,var,Q]
			]
		]


TakeDerivativeOfOdd[Letter_,SqrtRepl_,vars_,x_,CurrentSqrt_]:=
	Module[{P,q,p},
		q=(CurrentSqrt)^2/.SqrtRepl;
		P=(Letter)//Exp//Numerator//#/.{CurrentSqrt->0}&;
		p=-(Letter//Exp//Numerator//(#-P)/(CurrentSqrt)&);
		(2 P q D[p,x]-2 q p D[P,x]+P p D[q,x])/(CurrentSqrt*(-P^2+q p^2))/.SqrtRepl//Together
		]


TakeDerivativeOfOddModSqrt[Letter_,SqrtRepl_,vars_,x_,CurrentSqrt_]:=
	Module[{P,q,p},
		q=(CurrentSqrt)^2/.SqrtRepl;
		P=(Letter)//Exp//Numerator//#/.{CurrentSqrt->0}&;
		p=-(Letter//Exp//Numerator//(#-P)/(CurrentSqrt)&);
		(2 P q D[p,x]-2 q p D[P,x]+P p D[q,x])/(-P^2+p^2 q)/.SqrtRepl//Together
		]


(*SortAlphabet[EvenAlpha_,vars_,head_:Global`\[Alpha]Lett]:=Module[{lett},SortBy[EvenAlpha,{(Exp[#[[2]]]//ResourceFunction["PolynomialDegree"][#,vars]&)&,#[[1,1]]&}]//Table[{head[j]->#[[j]][[1]],head[j]->#[[j]][[2]]},{j,1,Length[#]}]&//Transpose];*)
(*SortAlphabet[EvenAlpha_,vars_,head_:Global`\[Alpha]Lett]:=Module[{lett},
SortBy[EvenAlpha,{Max[Cases[CoefficientRules[Exp[#[[2]]]],v_?VectorQ:>Total[v],2]]&,#[[1,1]]&}]//Table[{head[j]->#[[j]][[1]],head[j]->#[[j]][[2]]},{j,1,Length[#]}]&//Transpose];*)


End[]
(*Protect@@Names["Effortless`*"];*)
EndPackage[];
