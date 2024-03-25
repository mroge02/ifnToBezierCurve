(* ifnToBezierCurve` package *)

BeginPackage["ifnToBezierCurve`"];

(*see "https://mathematica.stackexchange.com/questions/270575/convert-interpolating-function-to-a-bezier-curve"*)

ifnToBezierCurve // ClearAll;

Begin["Private`"];


(*Bezier pts from derivatives at x0,x1*)hermiteToBezierPoints // ClearAll;
hermiteToBezierPoints[{{x0_}, f0__}, {{x1_}, f1__}] := 
  hermiteToBezierPoints[{x0, x1}, {f0}, {f1}];
hermiteToBezierPoints[{x0_, x1_}, f0_?VectorQ, f1_?VectorQ] := 
  Module[{dx, order0, rng0, order1, rng1, n, yy0, yy1}, dx = x1 - x0;
   order0 = Length@f0 - 1;
   rng0 = Range[0, order0];
   order1 = Length@f1 - 1;
   rng1 = Range[0, order1];
   n = order0 + order1 + 1;
   (*Bezier/Hermite order*)
   yy0 = f0*dx^rng0/Pochhammer[n + 1 - rng0, rng0];
   yy1 = f1*(-dx)^rng1/Pochhammer[n + 1 - rng1, rng1];
   Transpose@{Subdivide[x0, x1, n], 
     Join[Outer[Binomial, rng0, rng0] . yy0, 
      Outer[Binomial, rng1, rng1] . yy1 // Reverse]}];

(*"
 * Single BezierCurve per segment alternative
"*)

(*"unneeded as yet:
 *   $bitFlagsPos=2;(*bit field positions-inferred,perhaps mistaken*)
 *   $extrapBit=0;(*warn about extrapolation*)
 *   $rectArrayBit=1;(*data (f,f',...) is rect. array (not ragged)*)
 *   $machPrecBit=2;(*data (f,f',...) is MachinePrecision*)
 *   $repeatedBit=4;  (*repeated abscissae are permitted*)
"*)

ifnGetData // ClearAll;
iIfnGetData // ClearAll;
(***internal function***)
iIfnGetData[
    if : Verbatim[InterpolatingFunction][_, type_, 
      coords_, {Developer`PackedArrayForm, split_, a_}, ___], Automatic] /; 
   if["InterpolationMethod"] === "Hermite"(*/;BitAnd[type[[$bitFlagsPos]],
  2^$rectArrayBit]>
  0*):= <|"Coordinates" -> coords, 
   "FunctionData" -> Internal`PartitionRagged[a, Differences@split]|>;
iIfnGetData[
    if : Verbatim[InterpolatingFunction][_, type_, coords_, 
      a : {__List}, ___], Automatic] /; 
   if["InterpolationMethod"] === "Hermite" := <|"Coordinates" -> coords, 
   "FunctionData" -> a|>;
(*non-Hermite*)
iIfnGetData[if_InterpolatingFunction, Automatic] /; 
   if["InterpolationMethod"] =!= "Hermite" := (Message[ifnToBezierCurve::meth,
     if["InterpolationMethod"], First@if["InterpolationOrder"]];
   iIfnGetData[if, First@if["InterpolationOrder"]]);
(*non-Automatic SplineDegree*)
iIfnGetData[if_InterpolatingFunction, 
   sdeg_Integer?Positive] :=(*SplineDegree deg will be effectively be odd*)
  With[{coords = if@"Coordinates", 
    a = Transpose@
      Through[NestList[Derivative[1], if, Floor[sdeg/2]]["ValuesOnGrid"]]}, <|
    "Coordinates" -> coords, "FunctionData" -> a|>];
(***UI function***)
ifnGetData[if_InterpolatingFunction, sdeg_] := 
  Module[{data}, data = iIfnGetData[if, sdeg];
   With[{coordsegs = Split[First@data@"Coordinates", Unequal]}, {coordsegs, 
     Internal`PartitionRagged[data@"FunctionData", Length /@ coordsegs]}]];

ifnToBezierCurve // ClearAll;
ifnToBezierCurve::meth = 
  "Currently converting interpolations using interpolation method `` are \
unimplemented; returning a BezierCurve based on Hermite method with \
SplineDegree -> ``.";
(*If SplineDegree is not Automatic then Hermite interpolation model is used*)
ifnToBezierCurve // Options = {SplineDegree -> Automatic};
ifnToBezierCurve[ifn_InterpolatingFunction, OptionsPattern[]] /; 
   Length@ifn@"Domain" == 1 && ifn@"OutputDimensions" === {} := 
  With[{data = ifnGetData[ifn, OptionValue@SplineDegree]}, 
   MapThread[
     Function[{t, x},(*time steps,function step data*)
      If[MatrixQ[x] && 
        2  Length@x[[1]] > 
         First@ifn@
           "InterpolationOrder",(*all components the same sufficient degree*)
       With[{sd = 2  Length@x[[1]] - 1},(*SplineDegree*)
        BezierCurve[(*composite BezierCurve*)
         Join[{{t[[1]], x[[1, 1]]}}, 
          Flatten[MapThread[
            Rest[hermiteToBezierPoints[##]] &, {Partition[t, 2, 1], Most[x], 
             Rest[x]}], 1]], SplineDegree -> sd]
        ],
       (*degrees of components vary or are less than InterpolatingOrder*)
       Module[{delayedRest, iorder, ncoeffs, i},
        delayedRest := (delayedRest = Rest; Identity);
        ncoeffs = Length /@ x;
        JoinedCurve@(*joined BezierCurve pieces*)
         Table[iorder = ncoeffs[[j]] + ncoeffs[[j + 1]] - 1;
          i = 0;
          While[iorder < First@ifn@"InterpolationOrder", ++i;
           With[{dj = 3/4 - 1/4  (-1)^i  (3 + 2  i)}, 
            If[1 <= j + dj <= Length@x, iorder += ncoeffs[[j + dj]]]]];
          
          With[{bpts = 
             hermiteToBezierPoints[{t[[j]], t[[j + 1]]}, x[[j]], 
              Join[x[[j + 1]], 
               Table[Derivative[k][ifn][t[[j + 1]]], {k, ncoeffs[[j + 1]], 
                 iorder - ncoeffs[[j]]}]
               ]
              ]}, 
           BezierCurve[delayedRest@bpts, SplineDegree -> Length[bpts] - 1]
           ], {j, Length@t - 1}]
        ]
       ] (*end If[]*)
      ],
     data] /; FreeQ[data, ifnGetData]
   ];

End[];

EndPackage[];
