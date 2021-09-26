Sat Sep 25 10:53:57 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat3.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.10160560859        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -9.3670E+00 -7.8015E+01 -2.9921E+01 -7.5850E+01  3.8015E+01  5.8939E+00  3.0214E+00  3.9278E+00  1.7391E+00  2.6817E+01
            -1.6856E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1706.97470773549        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0230E+00  1.0860E+00  1.0333E+00  1.0183E+00  9.9248E-01  9.9860E-01  9.1725E-01  1.0475E+00  1.0066E+00  6.9527E-01
             1.0596E+00
 PARAMETER:  1.2274E-01  1.8252E-01  1.3275E-01  1.1811E-01  9.2450E-02  9.8596E-02  1.3628E-02  1.4638E-01  1.0659E-01 -2.6345E-01
             1.5793E-01
 GRADIENT:   4.2025E+01  1.1980E+01  7.2830E+00  2.2620E+01  1.3217E+01  4.9292E+00 -2.8685E+00 -8.7058E+00 -6.3314E-01 -8.2363E+00
            -3.5131E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1707.51912254440        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0259E+00  1.1873E+00  9.6566E-01  9.5844E-01  9.8955E-01  1.0062E+00  9.8846E-01  1.2737E+00  1.0105E+00  6.3975E-01
             1.0725E+00
 PARAMETER:  1.2555E-01  2.7168E-01  6.5061E-02  5.7547E-02  8.9490E-02  1.0620E-01  8.8391E-02  3.4192E-01  1.1044E-01 -3.4668E-01
             1.6996E-01
 GRADIENT:   4.6264E+01  3.3240E+01  5.3183E+00  2.4195E+01 -9.1997E+00  7.7162E+00  7.4016E+00  3.2005E-01 -1.9554E+00 -3.0362E+00
             4.0060E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.08277755300        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0147E+00  1.1799E+00  1.0175E+00  9.5610E-01  1.0156E+00  9.9703E-01  9.0074E-01  1.3207E+00  1.0426E+00  6.9177E-01
             1.0637E+00
 PARAMETER:  1.1455E-01  2.6540E-01  1.1734E-01  5.5112E-02  1.1550E-01  9.7029E-02 -4.5375E-03  3.7815E-01  1.4176E-01 -2.6850E-01
             1.6178E-01
 GRADIENT:   1.9417E+01  1.7382E+01  2.7690E+00  1.4569E+01 -2.6207E+00  3.7962E+00  1.5388E+00  1.8919E-01  2.7284E-01 -2.0766E+00
             1.6142E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1709.84545069016        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  1.0259E+00  1.5535E+00  6.7967E-01  7.1463E-01  1.0468E+00  1.0006E+00  7.9745E-01  1.0841E+00  1.2690E+00  7.1195E-01
             1.0579E+00
 PARAMETER:  1.2558E-01  5.4052E-01 -2.8614E-01 -2.3599E-01  1.4579E-01  1.0056E-01 -1.2634E-01  1.8074E-01  3.3824E-01 -2.3975E-01
             1.5626E-01
 GRADIENT:  -6.3543E+00  3.3703E+01  6.1233E+00  1.6329E+01 -1.4890E+01 -1.8592E+00  4.3577E+00 -3.6106E-01 -8.0841E-01  3.0199E-01
            -2.7870E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1711.81068083508        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      595
 NPARAMETR:  1.0319E+00  1.8136E+00  4.2362E-01  5.1859E-01  1.0776E+00  1.0114E+00  7.0671E-01  7.8546E-01  1.5264E+00  6.9922E-01
             1.0580E+00
 PARAMETER:  1.3144E-01  6.9531E-01 -7.5891E-01 -5.5664E-01  1.7477E-01  1.1129E-01 -2.4713E-01 -1.4149E-01  5.2291E-01 -2.5779E-01
             1.5641E-01
 GRADIENT:   5.9093E+00  8.2839E+00  1.7154E-01  5.6254E+00 -3.8665E+00  2.0743E+00 -7.8850E-01  5.1796E-01 -1.0536E+00  3.4235E-03
            -1.8226E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1711.88168403070        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.0290E+00  1.8163E+00  4.2119E-01  5.1194E-01  1.0811E+00  1.0058E+00  7.0893E-01  7.8224E-01  1.5301E+00  6.9880E-01
             1.0580E+00
 PARAMETER:  1.2858E-01  6.9683E-01 -7.6467E-01 -5.6954E-01  1.7794E-01  1.0579E-01 -2.4400E-01 -1.4559E-01  5.2532E-01 -2.5839E-01
             1.5634E-01
 GRADIENT:  -1.2360E-01 -2.8459E-02  4.1468E-01  1.4451E+00 -2.1318E+00 -2.0462E-02  3.0142E-02  4.9914E-01 -1.8474E+00 -1.5025E-01
            -1.9921E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1711.90258277545        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  1.0291E+00  1.8276E+00  4.2119E-01  5.0463E-01  1.0908E+00  1.0058E+00  7.0656E-01  7.8224E-01  1.5301E+00  7.0973E-01
             1.0590E+00
 PARAMETER:  1.2869E-01  7.0299E-01 -7.6467E-01 -5.8393E-01  1.8695E-01  1.0582E-01 -2.4735E-01 -1.4559E-01  5.2532E-01 -2.4287E-01
             1.5732E-01
 GRADIENT:   3.7188E-02  1.2755E-02  7.8315E-01 -1.5758E-02 -2.8875E-02  8.2978E-03 -1.5353E-02  3.9762E-01 -4.1529E+00 -9.2408E-03
             1.4815E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1711.90466903996        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  1.0291E+00  1.8275E+00  4.2109E-01  5.0455E-01  1.0909E+00  1.0058E+00  7.0661E-01  7.7972E-01  1.5302E+00  7.0990E-01
             1.0590E+00
 PARAMETER:  1.2869E-01  7.0299E-01 -7.6470E-01 -5.8393E-01  1.8695E-01  1.0582E-01 -2.4735E-01 -1.4886E-01  5.2551E-01 -2.4287E-01
             1.5732E-01
 GRADIENT:   8.0614E-03  1.6213E-01  1.2661E+05  3.3156E+05 -3.7795E-02  4.3926E-03 -1.3706E-02 -1.3007E+06  3.6840E+05 -2.0921E-02
            -7.0371E-03
 NUMSIGDIG:         4.5         4.1         3.3         3.3         4.0         4.0         3.2         3.3         3.3         2.7
                    4.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1102
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5010E-04 -2.8023E-02 -1.6010E-02  2.6708E-02 -3.8531E-02
 SE:             2.9861E-02  2.4504E-02  7.1305E-03  2.3480E-02  1.9258E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8797E-01  2.5280E-01  2.4749E-02  2.5534E-01  4.5420E-02

 ETASHRINKSD(%)  1.0000E-10  1.7907E+01  7.6112E+01  2.1339E+01  3.5482E+01
 ETASHRINKVR(%)  1.0000E-10  3.2607E+01  9.4294E+01  3.8124E+01  5.8374E+01
 EBVSHRINKSD(%)  4.6231E-01  1.8056E+01  7.7652E+01  2.2653E+01  3.5365E+01
 EBVSHRINKVR(%)  9.2249E-01  3.2852E+01  9.5005E+01  4.0174E+01  5.8223E+01
 RELATIVEINF(%)  9.9020E+01  4.8136E+00  3.9529E-01  4.2397E+00  8.1410E+00
 EPSSHRINKSD(%)  4.3831E+01
 EPSSHRINKVR(%)  6.8451E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1711.9046690399580     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -976.75384247621980     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1711.905       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.83E+00  4.21E-01  5.05E-01  1.09E+00  1.01E+00  7.07E-01  7.80E-01  1.53E+00  7.10E-01  1.06E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.03E+03
 
 TH 2
+       -7.22E+00  4.66E+02
 
 TH 3
+        3.98E+03 -6.11E+03  4.67E+08
 
 TH 4
+        4.59E+03 -6.70E+03 -2.80E+04  5.57E+08
 
 TH 5
+       -6.14E+00 -2.85E+02  3.91E+03  4.88E+03  7.90E+02
 
 TH 6
+       -2.82E+00 -1.14E+00  7.86E+03  9.06E+03 -9.06E-01  1.96E+02
 
 TH 7
+        4.85E-01  1.11E+01  2.54E+03  2.74E+03 -9.73E+00  8.00E-02  1.94E+02
 
 TH 8
+       -1.15E+04  1.77E+04 -9.78E+04 -1.42E+09 -1.16E+04 -2.27E+04 -7.01E+03  3.59E+09
 
 TH 9
+        1.52E+03 -2.54E+03  4.08E+04 -1.02E+04  1.85E+03  3.01E+03  1.05E+03 -4.16E+04  7.48E+07
 
 TH10
+        4.45E-01 -1.18E+01  1.97E+03  2.15E+03 -8.58E+01 -2.96E-01  2.50E+01 -5.50E+03  8.07E+02  8.90E+01
 
 TH11
+       -7.26E+00 -1.66E+01 -4.63E+05 -5.06E+05 -2.05E+01  1.75E+00  1.21E+01  1.28E+06 -1.85E+05  2.33E+01  1.90E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       20.014
Stop Time:
Sat Sep 25 10:54:18 CDT 2021
