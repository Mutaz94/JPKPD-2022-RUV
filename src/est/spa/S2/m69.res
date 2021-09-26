Sat Sep 25 12:26:20 CDT 2021
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
$DATA ../../../../data/spa/S2/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1659.48055427081        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.3673E+00  9.4682E+00 -7.7894E+01  1.2585E+02  1.2179E+02  4.9519E+00 -9.6673E+00  7.9713E+00 -2.4676E+00 -4.5156E+00
            -4.8779E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1673.80993868733        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.0253E+00  1.0569E+00  1.1326E+00  9.1091E-01  9.8998E-01  9.8467E-01  1.0647E+00  9.5078E-01  9.9217E-01  9.4629E-01
             1.1084E+00
 PARAMETER:  1.2495E-01  1.5538E-01  2.2455E-01  6.6912E-03  8.9934E-02  8.4547E-02  1.6269E-01  4.9525E-02  9.2135E-02  4.4793E-02
             2.0292E-01
 GRADIENT:   4.9364E+01  6.6689E+00  5.3681E-01  9.3715E+00  3.5092E+00  1.7041E-01 -2.7658E+00  3.7991E-01 -4.2430E-01 -4.2569E+00
            -8.5291E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1674.00657707785        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.0205E+00  1.0434E+00  1.0222E+00  9.0789E-01  9.4514E-01  9.9382E-01  1.1537E+00  7.6651E-01  9.6924E-01  9.2343E-01
             1.1098E+00
 PARAMETER:  1.2028E-01  1.4251E-01  1.2191E-01  3.3651E-03  4.3577E-02  9.3801E-02  2.4301E-01 -1.6591E-01  6.8761E-02  2.0340E-02
             2.0420E-01
 GRADIENT:   3.5860E+01  4.5315E+00  3.7012E+00 -2.0165E+00 -6.5336E+00  3.6047E+00  2.7406E+00  6.3505E-01  6.1083E-01 -2.7350E-01
             1.0048E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1674.30406412862        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  1.0309E+00  1.0953E+00  8.4726E-01  8.7154E-01  9.0335E-01  9.9637E-01  1.1148E+00  4.5832E-01  9.8876E-01  8.8810E-01
             1.1053E+00
 PARAMETER:  1.3048E-01  1.9102E-01 -6.5745E-02 -3.7493E-02 -1.6451E-03  9.6359E-02  2.0864E-01 -6.8018E-01  8.8700E-02 -1.8672E-02
             2.0009E-01
 GRADIENT:   6.0196E+00 -2.8747E-01 -3.1975E+00  2.6327E+00  1.0254E+00  1.1673E-01 -7.4481E-02  6.7626E-01  3.6756E-01  1.3908E+00
            -2.2658E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.50448519376        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.0296E+00  1.2185E+00  7.5620E-01  7.8935E-01  9.2144E-01  9.9607E-01  1.0291E+00  1.5237E-01  1.0589E+00  8.8584E-01
             1.1066E+00
 PARAMETER:  1.2913E-01  2.9763E-01 -1.7945E-01 -1.3655E-01  1.8186E-02  9.6061E-02  1.2866E-01 -1.7814E+00  1.5724E-01 -2.1219E-02
             2.0133E-01
 GRADIENT:   1.3699E+00 -1.3986E+00 -2.4196E+00  1.6307E+00  2.4730E+00 -1.6180E-01  3.4158E-02  6.1745E-02  1.4442E-01  4.7220E-01
             1.2764E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1674.73204811721        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      679
 NPARAMETR:  1.0290E+00  1.4752E+00  6.7905E-01  6.2793E-01  1.0141E+00  9.9553E-01  8.8180E-01  1.1552E-02  1.2704E+00  9.4794E-01
             1.1089E+00
 PARAMETER:  1.2862E-01  4.8880E-01 -2.8706E-01 -3.6533E-01  1.1398E-01  9.5515E-02 -2.5786E-02 -4.3609E+00  3.3931E-01  4.6540E-02
             2.0338E-01
 GRADIENT:  -2.7177E-01  3.4675E+00  4.9998E-01  1.5502E+00 -5.7230E+00 -2.3732E-01  7.1219E-02  2.9852E-04  5.2608E-01  5.8458E-01
             4.8758E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1675.01502761265        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      855
 NPARAMETR:  1.0287E+00  1.7379E+00  5.8218E-01  4.6028E-01  1.1263E+00  9.9580E-01  7.8343E-01  1.0000E-02  1.5848E+00  1.0227E+00
             1.1101E+00
 PARAMETER:  1.2831E-01  6.5268E-01 -4.4097E-01 -6.7593E-01  2.1894E-01  9.5789E-02 -1.4407E-01 -8.7721E+00  5.6046E-01  1.2243E-01
             2.0441E-01
 GRADIENT:  -1.2357E+00  8.0147E+00 -6.1609E-02  3.9756E+00  5.7721E-01 -2.2178E-02 -1.1537E+00  0.0000E+00 -3.2180E-01 -2.1878E-01
            -2.1489E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1675.07491348158        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1030
 NPARAMETR:  1.0289E+00  1.8266E+00  5.3430E-01  3.9721E-01  1.1621E+00  9.9599E-01  7.6515E-01  1.0000E-02  1.7379E+00  1.0460E+00
             1.1111E+00
 PARAMETER:  1.2849E-01  7.0245E-01 -5.2679E-01 -8.2328E-01  2.5025E-01  9.5981E-02 -1.6768E-01 -1.1094E+01  6.5269E-01  1.4502E-01
             2.0538E-01
 GRADIENT:  -7.4829E-01  7.3118E-01 -2.5345E-02  3.2328E-01  3.3633E-01  1.0674E-01  8.4801E-03  0.0000E+00 -4.3509E-02 -1.2389E-01
            -9.8380E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1675.07513139686        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  1.0292E+00  1.8321E+00  5.3105E-01  3.9319E-01  1.1641E+00  9.9570E-01  7.6365E-01  1.0000E-02  1.7491E+00  1.0479E+00
             1.1114E+00
 PARAMETER:  1.2883E-01  7.0544E-01 -5.3290E-01 -8.3347E-01  2.5193E-01  9.5689E-02 -1.6964E-01 -1.1255E+01  6.5908E-01  1.4681E-01
             2.0562E-01
 GRADIENT:   4.1519E-03  5.1957E-02 -4.2461E-03  3.1004E-02 -6.5561E-03 -6.7238E-04 -1.2368E-02  0.0000E+00  3.8723E-04 -6.6061E-03
            -5.9087E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1199
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2522E-04 -3.2847E-02 -1.8879E-04  3.1464E-02 -3.6981E-02
 SE:             2.9805E-02  2.3654E-02  6.8480E-05  2.1164E-02  2.3095E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9665E-01  1.6494E-01  5.8354E-03  1.3711E-01  1.0932E-01

 ETASHRINKSD(%)  1.4882E-01  2.0756E+01  9.9771E+01  2.9098E+01  2.2629E+01
 ETASHRINKVR(%)  2.9742E-01  3.7203E+01  9.9999E+01  4.9729E+01  4.0137E+01
 EBVSHRINKSD(%)  5.3728E-01  1.9145E+01  9.9802E+01  3.3294E+01  2.0092E+01
 EBVSHRINKVR(%)  1.0717E+00  3.4624E+01  1.0000E+02  5.5503E+01  3.6148E+01
 RELATIVEINF(%)  9.8856E+01  5.5044E+00  7.1985E-05  3.5939E+00  1.9806E+01
 EPSSHRINKSD(%)  4.3372E+01
 EPSSHRINKVR(%)  6.7933E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1675.0751313968599     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.92430483312171     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1675.075       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.83E+00  5.31E-01  3.93E-01  1.16E+00  9.96E-01  7.64E-01  1.00E-02  1.75E+00  1.05E+00  1.11E+00
 


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
+        1.05E+03
 
 TH 2
+       -8.36E+00  3.93E+02
 
 TH 3
+        5.80E+00  8.50E+01  1.96E+02
 
 TH 4
+       -2.02E+01  3.98E+02 -2.10E+02  1.08E+03
 
 TH 5
+       -6.80E+00 -1.29E+02 -1.78E+02  2.24E+02  4.79E+02
 
 TH 6
+       -6.03E-01 -1.40E+00  6.48E-01 -5.00E+00 -1.42E+00  1.94E+02
 
 TH 7
+        3.28E-01  2.08E+00 -3.32E+00 -1.88E+01 -1.49E+01  9.53E-01  1.66E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+00 -1.43E+01 -2.30E+01  5.92E+01  1.73E+00  9.92E-03  1.36E+01  0.00E+00  2.34E+01
 
 TH10
+       -3.18E+00 -1.51E+01 -2.36E+01  4.91E+00 -5.52E+01  3.75E+00  5.78E+00  0.00E+00  3.66E+00  8.08E+01
 
 TH11
+       -7.81E+00 -1.50E+01 -1.30E+01  4.32E+00  2.48E+00  2.37E+00  9.25E+00  0.00E+00  2.78E+00  1.60E+01  1.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.156
Stop Time:
Sat Sep 25 12:26:42 CDT 2021
