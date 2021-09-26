Sat Sep 25 07:57:38 CDT 2021
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
$DATA ../../../../data/spa/A1/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1287.72107600999        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.3508E+01 -2.7266E+01 -1.5033E+01 -5.1861E+01  3.4814E+01  3.3341E+01 -4.3059E+01  1.3154E+01 -7.5617E+01 -3.6457E+01
            -6.4325E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1470.41246042792        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0409E+00  9.2185E-01  1.2469E+00  1.0346E+00  1.0858E+00  9.2393E-01  1.3173E+00  8.0415E-01  1.2694E+00  1.1876E+00
             1.8591E+00
 PARAMETER:  1.4012E-01  1.8631E-02  3.2063E-01  1.3404E-01  1.8231E-01  2.0884E-02  3.7562E-01 -1.1798E-01  3.3852E-01  2.7190E-01
             7.2011E-01
 GRADIENT:   1.3371E+02 -3.9119E+01 -5.2395E+00 -6.3301E+01  1.1413E+01  4.4234E+00  2.0706E+00  4.3917E+00  1.5511E+01 -3.8142E+00
            -1.0081E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1475.32927120071        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0308E+00  6.4931E-01  1.2783E+00  1.3307E+00  1.0033E+00  8.7024E-01  1.7474E+00  4.7114E-01  9.7879E-01  1.2184E+00
             1.8607E+00
 PARAMETER:  1.3032E-01 -3.3185E-01  3.4557E-01  3.8573E-01  1.0334E-01 -3.8987E-02  6.5813E-01 -6.5261E-01  7.8558E-02  2.9758E-01
             7.2096E-01
 GRADIENT:   1.1842E+02  2.1011E+01 -1.6023E+01  6.9106E+01  2.1209E+01 -1.6682E+01  3.5988E+00  1.1623E+00 -1.5495E+01  5.2679E+00
            -1.9254E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1481.84025522800        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7879E-01  5.0489E-01  1.0361E+00  1.3558E+00  7.9226E-01  9.1205E-01  1.6321E+00  1.1945E-01  1.0222E+00  9.8260E-01
             1.9411E+00
 PARAMETER:  7.8564E-02 -5.8342E-01  1.3551E-01  4.0443E-01 -1.3286E-01  7.9378E-03  5.8988E-01 -2.0248E+00  1.2198E-01  8.2445E-02
             7.6326E-01
 GRADIENT:  -2.7921E+01  1.0572E+01  6.2115E+00  2.1055E+01 -1.1014E+01  4.9051E+00 -1.6829E+00  1.7895E-01  4.0008E-01  5.4613E-01
             7.4050E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1485.37550022689        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.8392E-01  1.6806E-01  7.8016E-01  1.4915E+00  5.7562E-01  8.9022E-01  2.8325E+00  1.0000E-02  9.4847E-01  8.0071E-01
             1.9213E+00
 PARAMETER:  8.3793E-02 -1.6834E+00 -1.4826E-01  4.9977E-01 -4.5232E-01 -1.6288E-02  1.1412E+00 -6.6661E+00  4.7092E-02 -1.2225E-01
             7.5300E-01
 GRADIENT:  -1.8755E+00  2.1034E+00  6.2172E-01  1.1678E+01 -4.3158E+00 -2.9172E+00 -3.4378E-01  0.0000E+00  1.5274E+00  6.3029E-01
            -5.8162E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1486.29661423990        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  9.8108E-01  5.6872E-02  7.7365E-01  1.5366E+00  5.5187E-01  8.9274E-01  4.6715E+00  1.0000E-02  9.1721E-01  7.8122E-01
             1.9597E+00
 PARAMETER:  8.0903E-02 -2.7670E+00 -1.5663E-01  5.2955E-01 -4.9444E-01 -1.3455E-02  1.6415E+00 -1.1369E+01  1.3577E-02 -1.4690E-01
             7.7279E-01
 GRADIENT:  -1.4367E+00  2.4045E-01  2.0475E+00 -2.1397E+00 -2.3201E+00 -6.8815E-01 -1.9788E-01  0.0000E+00 -5.5934E-01  6.2591E-02
             1.4378E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1486.59968282864        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.8014E-01  1.5050E-02  7.6545E-01  1.5575E+00  5.4121E-01  8.9457E-01  8.8750E+00  1.0000E-02  9.1008E-01  7.7580E-01
             1.9565E+00
 PARAMETER:  7.9939E-02 -4.0964E+00 -1.6730E-01  5.4307E-01 -5.1395E-01 -1.1414E-02  2.2832E+00 -1.7311E+01  5.7765E-03 -1.5386E-01
             7.7116E-01
 GRADIENT:  -2.8401E-01  5.4115E-02 -4.0835E-01  1.6781E+00  3.7445E-01  3.0009E-01 -4.9497E-02  0.0000E+00 -1.6123E-02  1.1024E-02
            -6.9420E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1487.14826767102        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      602
 NPARAMETR:  9.8386E-01  1.0000E-02  8.9621E-01  1.6039E+00  6.0758E-01  8.9467E-01  1.1898E+01  1.0000E-02  9.0475E-01  8.4825E-01
             1.9531E+00
 PARAMETER:  8.3727E-02 -4.6690E+00 -9.5838E-03  5.7245E-01 -3.9828E-01 -1.1303E-02  2.5763E+00 -1.9462E+01 -9.1414E-05 -6.4577E-02
             7.6944E-01
 GRADIENT:   5.7168E-01  0.0000E+00  5.0196E-01  1.9717E+00 -1.0474E+00 -1.2303E-01 -3.6801E-03  0.0000E+00  3.8074E-02  6.2189E-02
             7.1909E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1487.15005345032        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  9.8358E-01  1.0000E-02  8.9760E-01  1.6029E+00  6.0862E-01  8.9498E-01  1.1950E+01  1.0000E-02  9.0471E-01  8.4945E-01
             1.9501E+00
 PARAMETER:  8.3449E-02 -4.6726E+00 -8.0290E-03  5.7181E-01 -3.9656E-01 -1.0958E-02  2.5807E+00 -1.9477E+01 -1.4432E-04 -6.3169E-02
             7.6788E-01
 GRADIENT:  -2.5372E-03  0.0000E+00  2.2265E-03 -4.7534E-03 -2.3431E-03 -5.9933E-04 -1.9489E-04  0.0000E+00 -1.2021E-03  5.0188E-04
            -1.1664E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      764
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0381E-04  1.9874E-04 -7.7787E-06 -8.4642E-03 -2.0426E-02
 SE:             2.9318E-02  1.7762E-03  1.7235E-04  2.8197E-02  2.1932E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9445E-01  9.1091E-01  9.6400E-01  7.6404E-01  3.5168E-01

 ETASHRINKSD(%)  1.7820E+00  9.4049E+01  9.9423E+01  5.5376E+00  2.6526E+01
 ETASHRINKVR(%)  3.5322E+00  9.9646E+01  9.9997E+01  1.0768E+01  4.6015E+01
 EBVSHRINKSD(%)  1.8076E+00  9.4374E+01  9.9366E+01  5.0427E+00  2.5384E+01
 EBVSHRINKVR(%)  3.5826E+00  9.9684E+01  9.9996E+01  9.8312E+00  4.4324E+01
 RELATIVEINF(%)  8.4097E+01  9.4742E-03  2.4082E-04  4.6780E+00  1.9019E+00
 EPSSHRINKSD(%)  3.5495E+01
 EPSSHRINKVR(%)  5.8391E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1487.1500534503234     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -751.99922688658523     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.69
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1487.150       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.00E-02  8.98E-01  1.60E+00  6.09E-01  8.95E-01  1.19E+01  1.00E-02  9.05E-01  8.49E-01  1.95E+00
 


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
+        1.41E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.21E+01  0.00E+00  5.77E+02
 
 TH 4
+       -2.53E+01  0.00E+00 -4.97E+01  4.78E+02
 
 TH 5
+        1.86E+01  0.00E+00 -1.07E+03 -8.01E+01  2.30E+03
 
 TH 6
+       -2.42E+00  0.00E+00  3.11E+00 -7.10E+00  2.53E+00  2.45E+02
 
 TH 7
+       -2.02E-02  0.00E+00  6.31E-02 -1.41E-02  7.85E-04  7.22E-02  2.90E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.71E+01  0.00E+00  1.55E+01 -5.65E+00 -5.34E+00  2.71E+00 -1.92E-02  0.00E+00  2.13E+02
 
 TH10
+       -1.62E+01  0.00E+00 -2.22E+01 -4.19E+00 -9.34E+01 -1.90E+01  4.97E-02  0.00E+00  5.86E+00  1.08E+02
 
 TH11
+       -1.46E+01  0.00E+00 -9.02E+00 -6.85E+00  7.97E+00  4.88E+00  3.41E-03  0.00E+00  8.36E+00  1.92E+01  7.26E+01
 
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
 #CPUT: Total CPU Time in Seconds,       14.261
Stop Time:
Sat Sep 25 07:57:54 CDT 2021
