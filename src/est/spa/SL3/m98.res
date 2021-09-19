Sat Sep 18 13:10:29 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat98.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1625.72431618638        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0000E+02 -5.8635E+01 -4.6098E+01  5.6423E+00  1.1682E+02  2.5946E+01 -7.0455E+00 -4.4640E+00  1.3805E+01 -2.2931E+01
            -2.3815E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1633.94977380254        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.6017E-01  1.0682E+00  9.8829E-01  9.6564E-01  9.4573E-01  9.2777E-01  1.0467E+00  1.0526E+00  9.1390E-01  1.0399E+00
             1.0416E+00
 PARAMETER:  5.9358E-02  1.6602E-01  8.8224E-02  6.5039E-02  4.4197E-02  2.5024E-02  1.4561E-01  1.5124E-01  9.9607E-03  1.3915E-01
             1.4076E-01
 GRADIENT:   9.4268E-01 -3.9134E+00  8.7324E-02 -3.1197E+00  2.4013E+00 -5.6756E-01 -3.1068E+00 -2.4407E+00 -2.6086E+00  1.1260E+00
            -3.8643E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1634.65373670961        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  9.7315E-01  1.1818E+00  9.6735E-01  9.0502E-01  9.7401E-01  9.3668E-01  1.0013E+00  1.2706E+00  9.6560E-01  1.0029E+00
             1.0489E+00
 PARAMETER:  7.2787E-02  2.6702E-01  6.6807E-02  2.0448E-04  7.3662E-02  3.4584E-02  1.0135E-01  3.3949E-01  6.4997E-02  1.0293E-01
             1.4771E-01
 GRADIENT:  -3.3801E+00  4.9100E+00  5.7773E-01  6.6287E+00 -3.2578E+00 -5.0155E-01 -6.6171E-01  3.8307E-01  1.8242E-01 -2.9525E+00
            -5.5631E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.84076794883        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.7446E-01  1.2733E+00  9.5958E-01  8.4629E-01  1.0150E+00  9.3746E-01  9.5788E-01  1.3074E+00  1.0073E+00  1.0643E+00
             1.0481E+00
 PARAMETER:  7.4127E-02  3.4164E-01  5.8738E-02 -6.6892E-02  1.1488E-01  3.5419E-02  5.6968E-02  3.6806E-01  1.0724E-01  1.6229E-01
             1.4700E-01
 GRADIENT:  -3.1520E-01  5.9565E+00  2.3023E+00  3.1861E+00 -5.5850E+00 -1.8996E-01  2.4447E-01 -3.4018E-01 -4.4992E-01  5.0179E-01
            -7.3663E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1635.01418466570        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.7584E-01  1.5158E+00  7.1603E-01  6.8322E-01  1.0386E+00  9.3919E-01  8.6032E-01  1.2218E+00  1.1462E+00  1.0547E+00
             1.0496E+00
 PARAMETER:  7.5545E-02  5.1594E-01 -2.3403E-01 -2.8093E-01  1.3784E-01  3.7260E-02 -5.0453E-02  3.0032E-01  2.3644E-01  1.5323E-01
             1.4839E-01
 GRADIENT:  -6.6697E-02  6.9397E+00  4.8318E-01  4.1571E+00 -1.6660E+00  6.0902E-02 -7.4576E-02 -4.3539E-02 -2.8445E-02 -3.6612E-01
            -6.8736E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1635.04549666334        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.7608E-01  1.5882E+00  6.4527E-01  6.3016E-01  1.0503E+00  9.3923E-01  8.3614E-01  1.1770E+00  1.1998E+00  1.0610E+00
             1.0496E+00
 PARAMETER:  7.5794E-02  5.6261E-01 -3.3809E-01 -3.6178E-01  1.4908E-01  3.7302E-02 -7.8954E-02  2.6294E-01  2.8217E-01  1.5922E-01
             1.4840E-01
 GRADIENT:   2.7712E-02  4.0602E-01  1.3170E-02  3.0175E-01  5.2262E-02 -1.8654E-03 -7.2252E-02 -4.7205E-02 -6.8850E-02 -5.5312E-02
            -1.3596E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1635.04615648026        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      957            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7611E-01  1.6069E+00  6.3093E-01  6.1736E-01  1.0545E+00  9.3926E-01  8.3048E-01  1.1787E+00  1.2148E+00  1.0633E+00
             1.0496E+00
 PARAMETER:  7.5821E-02  5.7429E-01 -3.6056E-01 -3.8230E-01  1.5305E-01  3.7337E-02 -8.5755E-02  2.6441E-01  2.9459E-01  1.6138E-01
             1.4841E-01
 GRADIENT:   3.7250E+01  4.4008E+01  2.7131E-01  8.7405E+00  1.0248E+00  4.0262E+00  4.7703E-01  1.9284E-02  7.7805E-01  1.3232E-01
             1.1426E-01

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1635.04615760489        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  9.7611E-01  1.6069E+00  6.3090E-01  6.1737E-01  1.0545E+00  9.3926E-01  8.3047E-01  1.1790E+00  1.2148E+00  1.0633E+00
             1.0496E+00
 PARAMETER:  7.5819E-02  5.7429E-01 -3.6060E-01 -3.8228E-01  1.5306E-01  3.7335E-02 -8.5764E-02  2.6467E-01  2.9461E-01  1.6140E-01
             1.4841E-01
 GRADIENT:   8.4511E-03 -1.7295E-02 -5.6452E-03 -3.8835E-03  2.5529E-03  9.0304E-04  4.5489E-03  6.9505E-04  5.2196E-04  1.5534E-04
             7.0364E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1050
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2561E-04 -2.5576E-02 -2.8227E-02  2.0950E-02 -3.8233E-02
 SE:             2.9830E-02  2.3153E-02  1.0507E-02  2.1797E-02  2.2112E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9397E-01  2.6930E-01  7.2236E-03  3.3649E-01  8.3804E-02

 ETASHRINKSD(%)  6.5847E-02  2.2435E+01  6.4799E+01  2.6977E+01  2.5921E+01
 ETASHRINKVR(%)  1.3165E-01  3.9836E+01  8.7609E+01  4.6676E+01  4.5122E+01
 EBVSHRINKSD(%)  5.1674E-01  2.2108E+01  6.8592E+01  2.9169E+01  2.2707E+01
 EBVSHRINKVR(%)  1.0308E+00  3.9328E+01  9.0135E+01  4.9830E+01  4.0258E+01
 RELATIVEINF(%)  9.8824E+01  2.2580E+00  6.3899E-01  1.7426E+00  1.3066E+01
 EPSSHRINKSD(%)  4.5176E+01
 EPSSHRINKVR(%)  6.9943E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.0461576048913     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.89533104115310     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.046       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.61E+00  6.31E-01  6.17E-01  1.05E+00  9.39E-01  8.30E-01  1.18E+00  1.21E+00  1.06E+00  1.05E+00
 


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
+        1.30E+03
 
 TH 2
+       -7.71E+00  3.98E+02
 
 TH 3
+        9.25E+00  9.66E+01  2.18E+02
 
 TH 4
+       -2.04E+01  4.05E+02 -2.18E+02  9.91E+02
 
 TH 5
+       -3.74E+00 -1.60E+02 -2.59E+02  2.28E+02  5.35E+02
 
 TH 6
+       -3.10E-01 -1.41E+00  1.20E+00 -3.76E+00 -3.96E-01  2.20E+02
 
 TH 7
+       -1.48E+00  9.11E+00 -5.07E+00 -1.10E+01 -8.23E+00  2.28E+00  1.12E+02
 
 TH 8
+        3.26E-01 -9.31E+00 -1.85E+01  1.10E+01  2.13E+00 -4.49E-02  2.90E+00  4.59E+00
 
 TH 9
+        1.38E+00 -1.41E+01 -1.98E+01  4.55E+01 -3.06E-01 -8.20E-01  2.22E+01  3.67E+00  4.52E+01
 
 TH10
+       -6.13E-01 -1.18E+01 -2.17E+01  1.22E+00 -6.24E+01 -1.53E+00  9.83E+00  6.01E+00  4.51E+00  7.06E+01
 
 TH11
+       -7.12E+00 -1.26E+01 -1.40E+01  8.90E-01 -4.00E+00  4.02E+00  1.03E+01  3.97E+00  5.72E+00  1.04E+01  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.287
Stop Time:
Sat Sep 18 13:10:50 CDT 2021
