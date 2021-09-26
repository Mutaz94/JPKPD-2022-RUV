Sat Sep 25 08:00:55 CDT 2021
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
$DATA ../../../../data/spa/A1/dat35.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -949.943506480166        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.0444E+01  9.1684E+00  9.8259E+00 -1.0318E+01  9.7943E+01  1.3328E+01 -3.0847E+01  3.2472E+00 -4.4030E+01 -7.5985E+01
            -1.3397E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1429.15343837944        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0735E+00  9.8224E-01  1.0386E+00  1.0566E+00  9.7475E-01  8.7367E-01  9.9139E-01  9.3597E-01  9.9518E-01  9.3753E-01
             2.7666E+00
 PARAMETER:  1.7095E-01  8.2082E-02  1.3783E-01  1.5501E-01  7.4422E-02 -3.5049E-02  9.1353E-02  3.3829E-02  9.5164E-02  3.5498E-02
             1.1176E+00
 GRADIENT:   1.6974E+01 -2.1022E+00 -1.3880E+01  5.2027E+00  1.5118E+01 -1.5486E+01  3.8126E+00  5.7183E+00  5.0710E+00  1.0940E+01
             5.1366E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1434.34867402516        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0736E+00  7.8400E-01  7.0495E-01  1.1698E+00  6.9346E-01  9.4843E-01  7.7051E-01  4.5226E-01  9.5746E-01  4.4282E-01
             2.9062E+00
 PARAMETER:  1.7103E-01 -1.4335E-01 -2.4962E-01  2.5682E-01 -2.6606E-01  4.7048E-02 -1.6070E-01 -6.9349E-01  5.6524E-02 -7.1459E-01
             1.1668E+00
 GRADIENT:   2.2868E+00  2.7052E+01  3.5291E-01  4.9623E+01 -1.5839E+01  1.2009E+01 -4.3314E+00  1.9231E+00  2.5438E+00  1.5467E+00
             1.0109E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1440.53103648570        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0631E+00  5.8606E-01  6.0169E-01  1.2179E+00  5.6576E-01  8.8851E-01  1.5538E+00  5.6345E-02  8.1849E-01  4.3180E-01
             2.6758E+00
 PARAMETER:  1.6116E-01 -4.3434E-01 -4.0801E-01  2.9715E-01 -4.6959E-01 -1.8211E-02  5.4069E-01 -2.7763E+00 -1.0030E-01 -7.3980E-01
             1.0843E+00
 GRADIENT:  -1.2593E+01  1.0399E+01  1.7558E+01  4.5368E+00 -1.4419E+01 -1.3036E+01 -2.2818E+00  1.9528E-02 -4.1590E-01 -1.3385E+00
            -9.4981E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1443.80875900703        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0717E+00  4.6022E-01  3.3980E-01  1.1645E+00  3.7095E-01  9.2378E-01  1.7374E+00  1.0000E-02  8.0380E-01  3.5154E-01
             2.5475E+00
 PARAMETER:  1.6921E-01 -6.7606E-01 -9.7941E-01  2.5228E-01 -8.9168E-01  2.0723E-02  6.5238E-01 -5.1774E+00 -1.1840E-01 -9.4543E-01
             1.0351E+00
 GRADIENT:  -1.1131E+00  3.7673E+00 -2.0010E+00  1.2256E+01  2.4176E+00 -5.6384E+00  1.7623E+00  0.0000E+00 -7.8652E-01 -6.4018E-01
            -4.4285E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1444.21630126403        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      479
 NPARAMETR:  1.0767E+00  4.3436E-01  4.0034E-01  1.2140E+00  4.0069E-01  9.4092E-01  1.8242E+00  1.0000E-02  7.9991E-01  4.3631E-01
             2.5824E+00
 PARAMETER:  1.7394E-01 -7.3388E-01 -8.1545E-01  2.9391E-01 -8.1456E-01  3.9100E-02  7.0113E-01 -5.0641E+00 -1.2326E-01 -7.2939E-01
             1.0487E+00
 GRADIENT:   2.4156E+00  6.3707E+00  1.0261E+01  8.8875E-01 -1.6813E+01  2.3782E+00 -6.8247E-01  0.0000E+00 -6.3395E-01  2.7378E-01
             5.0412E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1445.01717505292        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  1.0690E+00  2.7303E-01  4.7515E-01  1.3218E+00  4.1800E-01  9.2929E-01  2.5318E+00  1.0000E-02  7.8913E-01  5.2630E-01
             2.5899E+00
 PARAMETER:  1.6671E-01 -1.1982E+00 -6.4412E-01  3.7897E-01 -7.7228E-01  2.6670E-02  1.0289E+00 -5.8265E+00 -1.3683E-01 -5.4188E-01
             1.0516E+00
 GRADIENT:   2.4535E-01  1.7163E+00  3.5877E+00  2.4377E-01 -4.9016E+00  6.7438E-01  9.6240E-01  0.0000E+00  1.6117E+00 -1.2619E+00
            -1.1716E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1445.41905910835        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0648E+00  1.8745E-01  4.8819E-01  1.3659E+00  4.1419E-01  9.2314E-01  3.1397E+00  1.0000E-02  7.6986E-01  5.7040E-01
             2.5853E+00
 PARAMETER:  1.6275E-01 -1.5742E+00 -6.1705E-01  4.1180E-01 -7.8142E-01  2.0030E-02  1.2441E+00 -6.7180E+00 -1.6155E-01 -4.6142E-01
             1.0498E+00
 GRADIENT:  -5.2102E-01  5.3848E-02  1.0852E-01  6.2432E+00 -1.4802E+00 -7.9230E-01 -1.1710E+00  0.0000E+00 -1.2978E+00  1.0649E+00
             1.4479E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1445.46529203924        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1006
 NPARAMETR:  1.0639E+00  1.6818E-01  4.7992E-01  1.3664E+00  4.0724E-01  9.2483E-01  3.3591E+00  1.0000E-02  7.7353E-01  5.5322E-01
             2.5877E+00
 PARAMETER:  1.6197E-01 -1.6827E+00 -6.3414E-01  4.1215E-01 -7.9835E-01  2.1859E-02  1.3117E+00 -7.0680E+00 -1.5680E-01 -4.9199E-01
             1.0508E+00
 GRADIENT:   2.1758E-02  1.3266E-02  2.0673E-02  2.7355E-02 -4.0926E-02  1.9436E-02  1.0584E-02  0.0000E+00  1.5876E-03 -2.2913E-03
            -3.0139E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1445.46529203924        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  1.0639E+00  1.6818E-01  4.7992E-01  1.3664E+00  4.0724E-01  9.2483E-01  3.3591E+00  1.0000E-02  7.7353E-01  5.5322E-01
             2.5877E+00
 PARAMETER:  1.6197E-01 -1.6827E+00 -6.3414E-01  4.1215E-01 -7.9835E-01  2.1859E-02  1.3117E+00 -7.0680E+00 -1.5680E-01 -4.9199E-01
             1.0508E+00
 GRADIENT:   2.1758E-02  1.3266E-02  2.0673E-02  2.7355E-02 -4.0926E-02  1.9436E-02  1.0584E-02  0.0000E+00  1.5876E-03 -2.2913E-03
            -3.0139E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1028
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5550E-03  3.9474E-02 -1.4777E-04 -2.6871E-02  6.7221E-03
 SE:             2.9068E-02  1.6278E-02  2.1111E-04  2.4735E-02  1.6585E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5734E-01  1.5308E-02  4.8393E-01  2.7733E-01  6.8524E-01

 ETASHRINKSD(%)  2.6191E+00  4.5466E+01  9.9293E+01  1.7133E+01  4.4439E+01
 ETASHRINKVR(%)  5.1695E+00  7.0261E+01  9.9995E+01  3.1331E+01  6.9130E+01
 EBVSHRINKSD(%)  2.7141E+00  5.9681E+01  9.9206E+01  1.4469E+01  4.0514E+01
 EBVSHRINKVR(%)  5.3545E+00  8.3744E+01  9.9994E+01  2.6844E+01  6.4614E+01
 RELATIVEINF(%)  9.2681E+01  6.3054E+00  1.9672E-04  2.1514E+01  1.1118E+00
 EPSSHRINKSD(%)  3.2124E+01
 EPSSHRINKVR(%)  5.3928E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1445.4652920392441     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -710.31446547550593     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1445.465       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.68E-01  4.80E-01  1.37E+00  4.07E-01  9.25E-01  3.36E+00  1.00E-02  7.74E-01  5.53E-01  2.59E+00
 


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
+        1.10E+03
 
 TH 2
+        4.00E+01  3.54E+03
 
 TH 3
+       -5.46E+01 -1.04E+03  3.84E+03
 
 TH 4
+       -8.21E+01 -3.28E+02 -3.04E+01  8.68E+02
 
 TH 5
+        1.79E+02  1.06E+03 -5.71E+03 -5.05E+02  9.23E+03
 
 TH 6
+       -4.81E+00  1.97E+01 -5.34E+00 -1.64E+01  2.85E+01  2.13E+02
 
 TH 7
+        1.05E+01  2.66E+02 -1.16E+02 -5.03E+01  1.59E+02  2.59E+00  2.19E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.87E+00 -2.74E+02  1.58E+02  3.43E+01 -1.73E+02 -6.90E+00 -2.05E+01  0.00E+00  2.32E+02
 
 TH10
+       -2.28E+01 -3.29E+02  5.21E+01  5.33E+01 -9.69E+01 -8.09E+00 -2.38E+01  0.00E+00  2.37E+01  1.19E+02
 
 TH11
+       -1.82E+01 -1.30E+02  4.40E+01  1.41E+01 -7.47E+01  3.02E+00 -9.44E+00  0.00E+00  2.16E+01  4.04E+01  4.98E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.133
Stop Time:
Sat Sep 25 08:01:17 CDT 2021
