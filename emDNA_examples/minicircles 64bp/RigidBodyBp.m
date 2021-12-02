(* ::Package:: *)

NumZero=10^-10;
XAxis={1,0,0};
YAxis={0,1,0};
ZAxis={0,0,1};


(* frame origin and matrix accessor *)
FrameAxesMatrix[frame_List]:=Return[frame[[2]]\[Transpose]];
FrameOrigin[frame_List]:=Return[frame[[1]]];
FrameMatrix[frame_List]:=Module[
	{mat},
	mat=DiagonalMatrix[{0,0,0,0}];
	mat[[1;;3,1;;3]]=FrameAxesMatrix[frame];
	mat[[1;;3,4]]=FrameOrigin[frame];
    mat[[4,4]]=1;
	Return[mat];
]


(* ZYZ Euler angles computation for two frames *)
EulerZYZAngles[frame1_List,frame2_List]:=Module[
	{Dmat,\[Kappa],\[Zeta],\[Eta]},

	(* rotation matrix *)
	Dmat=FrameAxesMatrix[frame1]\[Transpose].FrameAxesMatrix[frame2];
	
	(* \[Kappa] angle *)
	\[Kappa]=ArcCos[Dmat[[3,3]]];
	
	(* degenerated case *)
	If[
		\[Kappa]*\[Kappa]<NumZero*NumZero,
		Return[{ArcTan[Dmat[[1,1]],Dmat[[2,1]]],0,0}]
	];

	(* regular case *)
	\[Zeta]=ArcTan[Dmat[[1,3]],Dmat[[2,3]]];
	\[Eta]=ArcTan[-Dmat[[3,1]],Dmat[[3,2]]];
	
	(* minimal angle determination *)
	If[\[Zeta]+\[Eta]>\[Pi],\[Zeta]-=2\[Pi],If[\[Zeta]+\[Eta]<-\[Pi],\[Zeta]+=2\[Pi]]];

	Return[{\[Zeta],\[Kappa],\[Eta]}];

]
(* ZYZ Euler angles computation for anuglar rigid body parameters *)
EulerZYZAngles[\[Theta]1_,\[Theta]2_,\[Theta]3_]:=Module[
	{\[Kappa],\[Zeta],\[Eta]},
	
	(* degenerated case *)
	If[(\[Theta]1 *\[Theta]1*Degree^2 <NumZero*NumZero)&&(\[Theta]2 *\[Theta]2*Degree^2 <NumZero*NumZero),Return[{\[Theta]3 Degree,0,0}]];

	(* regular case *)
	\[Kappa]=Sqrt[(\[Theta]1*\[Theta]1+\[Theta]2*\[Theta]2)Degree*Degree];
	\[Zeta]=(\[Theta]3 Degree/2)-ArcTan[\[Theta]2 Degree,\[Theta]1 Degree];
	\[Eta]=(\[Theta]3 Degree/2)+ArcTan[\[Theta]2 Degree,\[Theta]1 Degree];

	Return[{\[Zeta],\[Kappa],\[Eta]}];

]


(* rotation matrix computation function *)
RotationMatrixD[\[Theta]1_,\[Theta]2_,\[Theta]3_]:=Module[
	{\[Zeta],\[Kappa],\[Eta]},
	{\[Zeta],\[Kappa],\[Eta]}=EulerZYZAngles@@(N[{\[Theta]1,\[Theta]2,\[Theta]3}]);
	Return[{{Cos[\[Zeta]] Cos[\[Eta]] Cos[\[Kappa]]-Sin[\[Zeta]] Sin[\[Eta]],-Cos[\[Eta]] Sin[\[Zeta]]-Cos[\[Zeta]] Cos[\[Kappa]] Sin[\[Eta]],Cos[\[Zeta]] Sin[\[Kappa]]},{Cos[\[Eta]] Cos[\[Kappa]] Sin[\[Zeta]]+Cos[\[Zeta]] Sin[\[Eta]],Cos[\[Zeta]] Cos[\[Eta]]-Cos[\[Kappa]] Sin[\[Zeta]] Sin[\[Eta]],Sin[\[Zeta]] Sin[\[Kappa]]},{-Cos[\[Eta]] Sin[\[Kappa]],Sin[\[Eta]] Sin[\[Kappa]],Cos[\[Kappa]]}}//N];
]


(* mid step rotation matrix computation function *)
MidStepRotationMatrix[\[Theta]1_,\[Theta]2_,\[Theta]3_]:=Module[
	{\[Zeta],\[Kappa],\[Eta]},
	{\[Zeta],\[Kappa],\[Eta]}=EulerZYZAngles@@(N[{\[Theta]1,\[Theta]2,\[Theta]3}]);
	Return[{{Cos[\[Zeta]] Cos[1/2 (-\[Zeta]+\[Eta])] Cos[\[Kappa]/2]-Sin[\[Zeta]] Sin[1/2 (-\[Zeta]+\[Eta])],-Cos[1/2 (-\[Zeta]+\[Eta])] Sin[\[Zeta]]-Cos[\[Zeta]] Cos[\[Kappa]/2] Sin[1/2 (-\[Zeta]+\[Eta])],Cos[\[Zeta]] Sin[\[Kappa]/2]},{Cos[1/2 (-\[Zeta]+\[Eta])] Cos[\[Kappa]/2] Sin[\[Zeta]]+Cos[\[Zeta]] Sin[1/2 (-\[Zeta]+\[Eta])],Cos[\[Zeta]] Cos[1/2 (-\[Zeta]+\[Eta])]-Cos[\[Kappa]/2] Sin[\[Zeta]] Sin[1/2 (-\[Zeta]+\[Eta])],Sin[\[Zeta]] Sin[\[Kappa]/2]},{-Cos[1/2 (-\[Zeta]+\[Eta])] Sin[\[Kappa]/2],Sin[1/2 (-\[Zeta]+\[Eta])] Sin[\[Kappa]/2],Cos[\[Kappa]/2]}}//N];
]


(* rigid body parameters functions *)
RigidBodyParameters[frame1_List,frame2_List]:=Module[
	{EulerAngles,\[Theta]1,\[Theta]2,\[Theta]3,midframe,\[Rho]1,\[Rho]2,\[Rho]3},
	
	(* Euler angles *)
	EulerAngles=EulerZYZAngles[frame1,frame2];

	(* angular parameters*)
	\[Theta]1=EulerAngles[[2]]Sin[(EulerAngles[[3]]-EulerAngles[[1]])/2]/Degree;
	\[Theta]2=EulerAngles[[2]]Cos[(EulerAngles[[3]]-EulerAngles[[1]])/2]/Degree;
	\[Theta]3=(EulerAngles[[3]]+EulerAngles[[1]])/Degree;

	(* displacement step parameters *)
	midframe=FrameAxesMatrix[frame1].MidStepRotationMatrix@@(N[{\[Theta]1,\[Theta]2,\[Theta]3}]);
	{\[Rho]1,\[Rho]2,\[Rho]3}=midframe\[Transpose].(FrameOrigin[frame2]-FrameOrigin[frame1]);

	(* rigid body parameters *)
	Return[{\[Theta]1,\[Theta]2,\[Theta]3,\[Rho]1,\[Rho]2,\[Rho]3}];

]
RigidBodyParametersList[frames_List]:=Module[
	{prms},
	prms=RigidBodyParameters@@#&/@Transpose[{Most[frames],Rest[frames]}];
	Return[prms];
]


(* frames rebuilding functions *)
RebuildFrame[RBprms_List,frame_List:{{0,0,0},{XAxis,YAxis,ZAxis}}]:=Module[
	{newframe,rotmat},

	newframe=frame;

	(* Euler angles *)
	rotmat=RotationMatrixD@@N[RBprms[[1;;3]]];

	(* new frame axis *)
	newframe[[2]]=(FrameAxesMatrix[frame].rotmat)\[Transpose];

	(* new frame origin *)
	rotmat=MidStepRotationMatrix@@N[RBprms[[1;;3]]];
	newframe[[1]]=(FrameAxesMatrix[frame].rotmat).RBprms[[4;;6]]+FrameOrigin[frame];

	Return[newframe];
	
]
RebuildFrames[RBprmsList_List,frame_List:{{0,0,0},{XAxis,YAxis,ZAxis}}]:=Module[
	{frames},
	frames=Table[frame,{i,Length[RBprmsList]+1}];

	(* loop *)
	Do[
		frames[[i+1]]=RebuildFrame[RBprmsList[[i]],frames[[i]]],
		{i,Length[RBprmsList]}
	];

	Return[frames];
	
]


(* frames transformation matrix function *)
FramesTransformationMatrix[frame1_List,frame2_List]:=Block[
	{rotmat,rvec},
	rotmat=FrameAxesMatrix[frame1]\[Transpose].FrameAxesMatrix[frame2];
	rvec=FrameAxesMatrix[frame1]\[Transpose].(FrameOrigin[frame2]-FrameOrigin[frame1]);
	Return[{rotmat,rvec}//Transpose//Flatten//Partition[#,4]&//Append[#,{0,0,0,1}]&];
];


(* framed tube rendering function *)
RenderFramedTube[framelist_List,closed_:False,lgaxis_Real:10.,shaxis_Real:3.,color_:LightBlue]:=Block[
	{bps,bpdimi,bpdimj,tube},
	bps={};
	bpdimi=shaxis;
	bpdimj=lgaxis;
	Do[
		AppendTo[
			bps,
			{
				framelist[[i]][[1]]-framelist[[i]][[2,1]]*bpdimi-framelist[[i]][[2,2]]*bpdimj,
				framelist[[i]][[1]]-framelist[[i]][[2,1]]*bpdimi+framelist[[i]][[2,2]]*bpdimj,
				framelist[[i]][[1]]+framelist[[i]][[2,1]]*bpdimi+framelist[[i]][[2,2]]*bpdimj,
				framelist[[i]][[1]]+framelist[[i]][[2,1]]*bpdimi-framelist[[i]][[2,2]]*bpdimj
			}
		],
		{i,1,Length[framelist]}
	];
	tube={
		Table[
			{bps[[i,1]],bps[[i+1,1]],bps[[i+1,2]],bps[[i,2]]},
			{i,1,Length[bps]-1}
		]//Map[Polygon,#]&,
		Table[
			{bps[[i,2]],bps[[i+1,2]],bps[[i+1,3]],bps[[i,3]]},
			{i,1,Length[bps]-1}
		]//Map[Polygon,#]&,
		Table[
			{bps[[i,3]],bps[[i+1,3]],bps[[i+1,4]],bps[[i,4]]},
			{i,1,Length[bps]-1}
		]//Map[Polygon,#]&,
		Table[
			{bps[[i,4]],bps[[i+1,4]],bps[[i+1,1]],bps[[i,1]]},
			{i,1,Length[bps]-1}
		]//Map[Polygon,#]&
	};
	If[
	Not[closed],
		AppendTo[tube,Polygon[{bps[[1,1]],bps[[1,2]],bps[[1,3]],bps[[1,4]]}]];
		AppendTo[
			tube,
			Polygon[{bps[[Length[bps],1]],bps[[Length[bps],2]],bps[[Length[bps],3]],bps[[Length[bps],4]]}]
		];
	];
	Return[tube];
];
RenderFrame[frame_List,lgaxis_Real:10.,shaxis_Real:3.]:=Module[
	{bps,poly,bpdimi,bpdimj},
	bpdimi=shaxis;
	bpdimj=lgaxis;
	bps={
		frame[[1]]-frame[[2,1]]*bpdimi-frame[[2,2]]*bpdimj,
		frame[[1]]-frame[[2,1]]*bpdimi+frame[[2,2]]*bpdimj,
		frame[[1]]+frame[[2,1]]*bpdimi+frame[[2,2]]*bpdimj,
		frame[[1]]+frame[[2,1]]*bpdimi-frame[[2,2]]*bpdimj
	};
	poly=Polygon[{bps[[1]],bps[[2]],bps[[3]],bps[[4]]}]
]
