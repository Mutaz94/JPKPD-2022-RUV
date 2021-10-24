Sat Oct 23 14:49:08 CDT 2021
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
$DATA ../../../../data/SD1/SL1/dat39.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3124.60075845701        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5820E+02  3.0723E+01  2.0307E+02  7.7923E+01 -1.3511E+00  5.4765E+01 -1.0245E+01 -1.4231E+02 -7.9310E+00  2.2291E+01
            -1.3108E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3386.43629817875        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.5062E-01  1.0723E+00  8.1319E-01  9.7114E-01  1.0116E+00  9.4177E-01  9.6412E-01  1.1254E+00  9.7165E-01  9.2689E-01
             1.6175E+00
 PARAMETER:  4.9357E-02  1.6981E-01 -1.0679E-01  7.0713E-02  1.1157E-01  4.0008E-02  6.3455E-02  2.1816E-01  7.1236E-02  2.4085E-02
             5.8089E-01
 GRADIENT:   5.8047E+01  2.4763E+01 -1.6004E+01  1.1354E+01  5.4749E-01  7.4164E+00  1.0191E+01 -3.0540E+00  4.6734E+00  3.9905E+00
             1.3133E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3389.05646714375        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.4498E-01  1.0840E+00  8.2138E-01  9.7107E-01  1.0374E+00  9.5320E-01  9.4635E-01  1.7213E+00  9.1208E-01  1.0558E+00
             1.5615E+00
 PARAMETER:  4.3414E-02  1.8069E-01 -9.6764E-02  7.0648E-02  1.3671E-01  5.2068E-02  4.4856E-02  6.4310E-01  7.9708E-03  1.5430E-01
             5.4562E-01
 GRADIENT:   5.8714E+01  1.8803E+01 -1.5498E+01  2.2359E+01  7.8354E+00  1.2611E+01  1.2014E+01  1.7402E+01 -5.2796E+00  9.3893E+00
             1.2032E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3410.32764224574        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  9.9693E-01  1.4350E+00  9.8773E-01  7.9583E-01  1.3701E+00  9.5998E-01  6.4442E-01  2.0352E+00  1.0830E+00  1.2593E+00
             1.4708E+00
 PARAMETER:  9.6922E-02  4.6114E-01  8.7653E-02 -1.2837E-01  4.1489E-01  5.9158E-02 -3.3940E-01  8.1060E-01  1.7973E-01  3.3053E-01
             4.8582E-01
 GRADIENT:   1.3432E+01  8.7250E+00  1.4835E+01  3.0920E+01 -6.0276E+00  7.6805E-01 -2.2673E+00 -2.6221E+00 -2.9948E+00  1.9286E-01
            -2.9197E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3417.53296261666        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  9.8872E-01  1.7983E+00  8.5983E-01  5.6580E-01  1.6974E+00  9.5575E-01  6.7710E-01  2.3160E+00  1.2260E+00  1.3461E+00
             1.5103E+00
 PARAMETER:  8.8660E-02  6.8686E-01 -5.1021E-02 -4.6951E-01  6.2912E-01  5.4746E-02 -2.8994E-01  9.3984E-01  3.0380E-01  3.9719E-01
             5.1228E-01
 GRADIENT:  -7.4689E+00  2.5902E+01  1.7748E+00  2.0365E+01  9.9141E+00 -1.2176E+00 -4.5010E+00 -3.0555E+00  4.6369E+00 -1.2282E+00
            -3.0654E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3417.74754029126        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  9.8768E-01  1.7892E+00  8.5966E-01  5.5859E-01  1.7069E+00  9.4785E-01  7.0656E-01  2.3129E+00  1.2208E+00  1.3451E+00
             1.5167E+00
 PARAMETER:  8.7608E-02  6.8179E-01 -5.1215E-02 -4.8234E-01  6.3469E-01  4.6437E-02 -2.4735E-01  9.3849E-01  2.9948E-01  3.9648E-01
             5.1651E-01
 GRADIENT:   1.8022E+02  5.2267E+02  2.7671E+00  6.6822E+01  1.0805E+02  1.4085E+01  9.4632E+00 -2.1141E-01  1.3000E+01  8.7546E+00
             1.8607E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3417.95771524120        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.8811E-01  1.7981E+00  8.5949E-01  5.5299E-01  1.6896E+00  9.5950E-01  7.1638E-01  2.3306E+00  1.2183E+00  1.3468E+00
             1.5103E+00
 PARAMETER:  8.8043E-02  6.8671E-01 -5.1415E-02 -4.9241E-01  6.2450E-01  5.8661E-02 -2.3354E-01  9.4613E-01  2.9742E-01  3.9774E-01
             5.1232E-01
 GRADIENT:  -8.4833E+00  2.0896E-01  4.3615E+00  3.4401E+00 -8.6290E-01  4.0721E-01  2.3303E+00 -2.7826E+00  8.7399E+00 -1.6966E+00
             4.0192E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3418.10086752884        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  9.8880E-01  1.8003E+00  8.5485E-01  5.5099E-01  1.6867E+00  9.5868E-01  7.0497E-01  2.3345E+00  1.2030E+00  1.3486E+00
             1.5091E+00
 PARAMETER:  8.8732E-02  6.8793E-01 -5.6832E-02 -4.9604E-01  6.2280E-01  5.7805E-02 -2.4960E-01  9.4779E-01  2.8479E-01  3.9903E-01
             5.1154E-01
 GRADIENT:  -6.8554E+00  5.4490E-01  4.2787E+00  1.8466E+00 -3.5237E+00  8.8105E-02 -1.1386E+00 -2.6641E+00  6.1780E+00 -1.6580E+00
            -1.6938E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -3418.16456292521        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  9.8876E-01  1.8005E+00  8.5488E-01  5.5088E-01  1.6871E+00  9.5853E-01  7.0558E-01  2.3337E+00  1.1897E+00  1.3483E+00
             1.5094E+00
 PARAMETER:  8.8693E-02  6.8809E-01 -5.6793E-02 -4.9625E-01  6.2303E-01  5.7650E-02 -2.4874E-01  9.4744E-01  2.7369E-01  3.9888E-01
             5.1173E-01
 GRADIENT:   5.1581E+02 -7.3580E+01 -5.1793E+02  1.0750E+02 -8.6690E+01 -1.9317E-02 -1.0715E+02  5.3863E+01 -1.8589E+02  1.2927E+02
            -1.0599E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1102
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1506E-03 -3.5533E-02 -3.2499E-02  3.7162E-02 -2.8665E-02
 SE:             2.9778E-02  2.3291E-02  1.7423E-02  2.1605E-02  2.6224E-02
 N:                     100         100         100         100         100

 P VAL.:         8.8915E-01  1.2710E-01  6.2138E-02  8.5420E-02  2.7436E-01

 ETASHRINKSD(%)  2.3885E-01  2.1974E+01  4.1632E+01  2.7620E+01  1.2146E+01
 ETASHRINKVR(%)  4.7713E-01  3.9119E+01  6.5932E+01  4.7612E+01  2.2817E+01
 EBVSHRINKSD(%)  6.0973E-01  2.1172E+01  4.5662E+01  3.1140E+01  1.0555E+01
 EBVSHRINKVR(%)  1.2157E+00  3.7861E+01  7.0474E+01  5.2583E+01  1.9995E+01
 RELATIVEINF(%)  9.8776E+01  8.5512E+00  1.6990E+01  6.2634E+00  2.9485E+01
 EPSSHRINKSD(%)  1.9743E+01
 EPSSHRINKVR(%)  3.5587E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3418.1645629252071     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1764.0752031567963     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3418.165       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.80E+00  8.55E-01  5.51E-01  1.69E+00  9.59E-01  7.06E-01  2.33E+00  1.19E+00  1.35E+00  1.51E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       80.755
Stop Time:
Sat Oct 23 14:49:22 CDT 2021
