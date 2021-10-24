Sat Oct 23 23:20:15 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD3/SL1/dat11.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m11.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2085.55138725643        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4594E+02 -5.9705E+00 -4.3565E+01  6.2411E+01  7.2301E+01  6.1263E+01 -3.6258E+00  9.1357E+00  3.2558E+00  1.0391E+01
             4.2955E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2090.02954710340        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9740E-01  1.0372E+00  1.0622E+00  1.0133E+00  9.7640E-01  9.4399E-01  1.0272E+00  9.6828E-01  1.0219E+00  9.5133E-01
             9.2333E-01
 PARAMETER:  9.7392E-02  1.3651E-01  1.6030E-01  1.1318E-01  7.6120E-02  4.2363E-02  1.2688E-01  6.7766E-02  1.2170E-01  5.0106E-02
             2.0228E-02
 GRADIENT:   1.7334E+02  2.5106E+01  5.2901E+00  2.8277E+01 -1.3048E+01 -1.7309E+00 -3.6759E+00  2.3822E+00 -1.3545E-01  2.6058E+00
            -2.4475E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2091.59015512069        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8545E-01  9.1193E-01  1.0158E+00  1.1051E+00  8.9488E-01  9.4220E-01  1.2760E+00  7.4600E-01  8.9952E-01  8.1186E-01
             9.4945E-01
 PARAMETER:  8.5347E-02  7.8125E-03  1.1569E-01  1.9995E-01 -1.1069E-02  4.0463E-02  3.4372E-01 -1.9303E-01 -5.8984E-03 -1.0843E-01
             4.8126E-02
 GRADIENT:   1.4528E+02  4.3264E+01  2.2198E+01  5.7526E+01 -2.0305E+01  1.4027E+00 -3.0955E+00 -4.8107E+00 -8.5633E+00 -1.1736E+01
            -4.6679E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2098.05475952890        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.2983E-01  9.3590E-01  8.9119E-01  1.0418E+00  8.6704E-01  9.1621E-01  1.2497E+00  6.2660E-01  9.3557E-01  8.3709E-01
             9.4866E-01
 PARAMETER:  2.7249E-02  3.3756E-02 -1.5201E-02  1.4096E-01 -4.2674E-02  1.2493E-02  3.2289E-01 -3.6744E-01  3.3400E-02 -7.7824E-02
             4.7300E-02
 GRADIENT:  -1.1076E+00  3.0234E+00  3.1197E+00  4.8675E-01 -4.6310E+00 -1.3242E+00  2.3373E-02 -4.9303E-01 -7.8313E-01 -8.9875E-02
            -1.4565E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2098.07530887326        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.2993E-01  8.9424E-01  9.4074E-01  1.0701E+00  8.7555E-01  9.1851E-01  1.2787E+00  6.9725E-01  9.2409E-01  8.5154E-01
             9.4850E-01
 PARAMETER:  2.7353E-02 -1.1786E-02  3.8914E-02  1.6777E-01 -3.2899E-02  1.4997E-02  3.4587E-01 -2.6061E-01  2.1057E-02 -6.0704E-02
             4.7127E-02
 GRADIENT:   6.6342E-01  5.0776E-01 -3.2990E-01  6.2498E-01 -2.7266E-01 -1.2304E-01 -2.7699E-01  1.3578E-01 -4.1608E-01  3.7739E-01
            -7.3584E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2098.07537847625        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.2978E-01  8.8155E-01  9.5440E-01  1.0790E+00  8.7657E-01  9.1844E-01  1.2897E+00  7.1187E-01  9.1991E-01  8.5455E-01
             9.4850E-01
 PARAMETER:  2.7192E-02 -2.6069E-02  5.3323E-02  1.7603E-01 -3.1742E-02  1.4919E-02  3.5444E-01 -2.3985E-01  1.6516E-02 -5.7176E-02
             4.7125E-02
 GRADIENT:   6.8278E-01  7.2545E-01 -2.1269E-01  8.4496E-01 -5.7530E-01 -9.1238E-02 -3.0022E-01  1.4990E-01 -4.4437E-01  4.1108E-01
            -7.9512E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2098.07541750880        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.2965E-01  8.7158E-01  9.6492E-01  1.0859E+00  8.7730E-01  9.1837E-01  1.2988E+00  7.2273E-01  9.1672E-01  8.5670E-01
             9.4850E-01
 PARAMETER:  2.7053E-02 -3.7453E-02  6.4286E-02  1.8240E-01 -3.0904E-02  1.4845E-02  3.6142E-01 -2.2472E-01  1.3042E-02 -5.4667E-02
             4.7131E-02
 GRADIENT:   6.6724E-01  8.3333E-01 -1.3604E-01  9.4881E-01 -7.3285E-01 -7.1575E-02 -3.0264E-01  1.5257E-01 -4.4598E-01  4.1574E-01
            -8.0489E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2098.07542776388        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.2956E-01  8.6440E-01  9.7239E-01  1.0908E+00  8.7781E-01  9.1832E-01  1.3054E+00  7.3028E-01  9.1446E-01  8.5817E-01
             9.4851E-01
 PARAMETER:  2.6951E-02 -4.5724E-02  7.2001E-02  1.8693E-01 -3.0329E-02  1.4786E-02  3.6652E-01 -2.1433E-01  1.0580E-02 -5.2954E-02
             4.7136E-02
 GRADIENT:   6.5053E-01  8.8117E-01 -9.2458E-02  9.9211E-01 -8.0803E-01 -5.9408E-02 -2.9971E-01  1.5183E-01 -4.4089E-01  4.1245E-01
            -7.9997E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2098.07856501395        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1374
 NPARAMETR:  9.2977E-01  8.6334E-01  9.7253E-01  1.0904E+00  8.7801E-01  9.1877E-01  1.3114E+00  7.2680E-01  9.1605E-01  8.5525E-01
             9.4865E-01
 PARAMETER:  2.7179E-02 -4.6945E-02  7.2147E-02  1.8652E-01 -3.0097E-02  1.5277E-02  3.7113E-01 -2.1911E-01  1.2317E-02 -5.6363E-02
             4.7285E-02
 GRADIENT:   1.2856E+00  6.3584E-02  1.2273E-01 -4.7865E-01  1.2450E-01  1.3619E-01  8.2540E-02 -4.2384E-02  1.4167E-01 -7.4504E-02
            -5.3261E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1374
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4942E-04 -1.3167E-03 -2.5489E-02 -1.9829E-03 -2.0746E-02
 SE:             2.9883E-02  1.9763E-02  1.2223E-02  2.5085E-02  2.2138E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8533E-01  9.4688E-01  3.7045E-02  9.3700E-01  3.4869E-01

 ETASHRINKSD(%)  1.0000E-10  3.3790E+01  5.9051E+01  1.5961E+01  2.5836E+01
 ETASHRINKVR(%)  1.0000E-10  5.6163E+01  8.3231E+01  2.9374E+01  4.4996E+01
 EBVSHRINKSD(%)  3.5542E-01  3.4375E+01  6.1482E+01  1.6014E+01  2.4401E+01
 EBVSHRINKVR(%)  7.0959E-01  5.6934E+01  8.5163E+01  2.9464E+01  4.2848E+01
 RELATIVEINF(%)  9.8584E+01  2.7436E+00  1.9734E+00  5.3056E+00  9.5474E+00
 EPSSHRINKSD(%)  3.4166E+01
 EPSSHRINKVR(%)  5.6659E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2098.0785650139492     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1179.1400318092765     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2098.079       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.30E-01  8.63E-01  9.73E-01  1.09E+00  8.78E-01  9.19E-01  1.31E+00  7.27E-01  9.16E-01  8.55E-01  9.49E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      109.781
Stop Time:
Sat Oct 23 23:20:32 CDT 2021
