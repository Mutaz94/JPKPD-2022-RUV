Sat Sep 18 09:59:14 CDT 2021
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
$DATA ../../../../data/spa/A2/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -690.065791514912        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7839E-01  3.5569E+01  2.3016E+01  5.3801E+01  2.7137E+02  5.0847E+01 -6.2227E+01 -4.0124E+01 -7.0600E+01 -1.7749E+02
            -1.6653E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1336.79526874586        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0239E+00  9.3386E-01  9.3974E-01  1.0282E+00  8.4979E-01  8.1667E-01  1.0532E+00  1.0016E+00  1.0796E+00  1.0686E+00
             2.2920E+00
 PARAMETER:  1.2359E-01  3.1572E-02  3.7848E-02  1.2778E-01 -6.2761E-02 -1.0252E-01  1.5179E-01  1.0163E-01  1.7660E-01  1.6634E-01
             9.2944E-01
 GRADIENT:  -1.9417E+01  3.7973E+00  1.0685E+00  2.9214E+00  5.9096E+01 -1.1149E+01 -1.0927E+01 -2.4102E+00 -1.3813E+01 -2.3764E+01
            -1.6753E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1348.71523899879        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0078E+00  5.2939E-01  2.8083E-01  1.2185E+00  3.2777E-01  8.3992E-01  1.4325E+00  4.0944E-01  1.1395E+00  4.2554E-01
             2.1924E+00
 PARAMETER:  1.0775E-01 -5.3603E-01 -1.1700E+00  2.9762E-01 -1.0154E+00 -7.4451E-02  4.5940E-01 -7.9296E-01  2.3062E-01 -7.5439E-01
             8.8498E-01
 GRADIENT:  -1.0051E+02  4.4817E+01 -7.9845E+00  2.1237E+02 -7.1812E+00 -1.7769E+01 -2.4453E+00 -8.2194E+00 -5.3146E+00 -2.1826E+01
            -1.3743E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1380.11427160133        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0479E+00  5.0720E-01  3.2491E-01  1.1369E+00  3.4593E-01  8.7414E-01  1.2565E+00  1.9502E-01  1.0116E+00  7.2022E-01
             2.6368E+00
 PARAMETER:  1.4680E-01 -5.7885E-01 -1.0242E+00  2.2828E-01 -9.6153E-01 -3.4514E-02  3.2831E-01 -1.5346E+00  1.1150E-01 -2.2820E-01
             1.0696E+00
 GRADIENT:   1.0295E-01  2.4476E+01  4.4314E+00  4.3725E+01 -9.5626E+00  5.4648E+00 -5.1427E-01  1.5264E-03 -6.5086E+00  7.2040E+00
            -1.1552E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1382.30506355916        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  1.0489E+00  3.8851E-01  2.6106E-01  1.1029E+00  2.8204E-01  8.5772E-01  1.2962E+00  7.1794E-02  1.0570E+00  6.7826E-01
             2.6879E+00
 PARAMETER:  1.4778E-01 -8.4543E-01 -1.2430E+00  1.9793E-01 -1.1657E+00 -5.3478E-02  3.5948E-01 -2.5339E+00  1.5542E-01 -2.8823E-01
             1.0888E+00
 GRADIENT:  -1.0854E+01  7.9659E+00  1.8465E+00  1.4446E+01 -2.0459E+01 -2.2530E+00  7.0903E-01 -3.6481E-02  2.2534E-01 -3.4717E-01
             5.6749E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1383.85029567112        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.0478E+00  2.7032E-01  3.5878E-01  1.2336E+00  3.1785E-01  8.5267E-01  1.7885E+00  4.2752E-02  9.8569E-01  7.6716E-01
             2.6857E+00
 PARAMETER:  1.4671E-01 -1.2081E+00 -9.2505E-01  3.0997E-01 -1.0462E+00 -5.9378E-02  6.8136E-01 -3.0523E+00  8.5586E-02 -1.6506E-01
             1.0879E+00
 GRADIENT:  -9.0657E-01  6.8497E+00  1.5773E+01  1.3590E+01 -2.8242E+01  1.0134E+00 -4.2549E-01 -3.9886E-04 -1.2620E+00  1.2727E+00
            -4.5345E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1385.39894071946        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  1.0364E+00  1.1350E-01  3.9290E-01  1.3167E+00  3.2362E-01  8.4394E-01  3.0619E+00  1.0000E-02  9.5224E-01  7.8895E-01
             2.7116E+00
 PARAMETER:  1.3579E-01 -2.0760E+00 -8.3420E-01  3.7514E-01 -1.0282E+00 -6.9669E-02  1.2190E+00 -7.0291E+00  5.1059E-02 -1.3705E-01
             1.0976E+00
 GRADIENT:  -5.9356E+00  5.4814E-01 -1.3303E-02  1.2737E+01 -2.5794E+00  7.0254E-01 -9.8026E-01  0.0000E+00  1.5354E-01  3.3141E-01
            -2.5602E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1385.70459841068        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      851
 NPARAMETR:  1.0339E+00  4.5610E-02  3.9677E-01  1.3375E+00  3.1926E-01  8.3882E-01  5.6756E+00  1.0000E-02  9.3837E-01  7.9274E-01
             2.7305E+00
 PARAMETER:  1.3334E-01 -2.9876E+00 -8.2441E-01  3.9081E-01 -1.0418E+00 -7.5764E-02  1.8362E+00 -1.1759E+01  3.6387E-02 -1.3226E-01
             1.1045E+00
 GRADIENT:  -9.9992E-01  4.5134E-01  4.7664E-01  1.3427E+00 -1.1495E+00 -1.5663E-01  6.8720E-01  0.0000E+00 -1.0715E-01 -1.7798E-01
             1.0403E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1385.73590039127        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  1.0326E+00  2.0771E-02  4.0229E-01  1.3502E+00  3.1988E-01  8.3861E-01  8.4018E+00  1.0000E-02  9.3332E-01  7.9809E-01
             2.7280E+00
 PARAMETER:  1.3207E-01 -3.7742E+00 -8.1058E-01  4.0025E-01 -1.0398E+00 -7.6011E-02  2.2285E+00 -1.6092E+01  3.0996E-02 -1.2553E-01
             1.1036E+00
 GRADIENT:   7.1915E-02  1.1566E-01  9.7300E-01  8.2507E-01 -1.5568E+00  6.7724E-02  1.6094E-01  0.0000E+00 -1.4219E-01 -2.9717E-02
            -5.2830E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1385.74626407923        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1204
 NPARAMETR:  1.0318E+00  1.0000E-02  4.0398E-01  1.3552E+00  3.2000E-01  8.3798E-01  1.2008E+01  1.0000E-02  9.3188E-01  7.9932E-01
             2.7317E+00
 PARAMETER:  1.3127E-01 -4.5138E+00 -8.0640E-01  4.0395E-01 -1.0394E+00 -7.6766E-02  2.5855E+00 -2.0272E+01  2.9444E-02 -1.2399E-01
             1.1049E+00
 GRADIENT:  -3.1271E-01  0.0000E+00  3.4259E-02  3.9645E-01 -1.4260E-01 -2.6880E-02  1.5587E-02  0.0000E+00  2.8863E-02  2.9764E-02
             1.2638E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1385.74638224822        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0319E+00  1.0000E-02  4.0358E-01  1.3546E+00  3.1982E-01  8.3806E-01  1.1972E+01  1.0000E-02  9.3191E-01  7.9913E-01
             2.7309E+00
 PARAMETER:  1.3138E-01 -4.5100E+00 -8.0739E-01  4.0348E-01 -1.0400E+00 -7.6662E-02  2.5826E+00 -2.0253E+01  2.9478E-02 -1.2424E-01
             1.1046E+00
 GRADIENT:   5.9817E-02  0.0000E+00 -4.9384E-02 -1.1315E-01  1.0032E-01  3.4484E-03  6.3552E-03  0.0000E+00  2.0912E-03 -7.4624E-03
            -1.2083E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1261
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.0950E-04  1.5850E-03  4.1120E-05 -9.5251E-03 -4.3548E-03
 SE:             2.8713E-02  1.8388E-03  2.1080E-04  2.6864E-02  2.2484E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8029E-01  3.8871E-01  8.4534E-01  7.2291E-01  8.4643E-01

 ETASHRINKSD(%)  3.8076E+00  9.3840E+01  9.9294E+01  1.0003E+01  2.4675E+01
 ETASHRINKVR(%)  7.4702E+00  9.9621E+01  9.9995E+01  1.9005E+01  4.3261E+01
 EBVSHRINKSD(%)  3.6728E+00  9.4630E+01  9.9306E+01  8.9682E+00  2.3849E+01
 EBVSHRINKVR(%)  7.2108E+00  9.9712E+01  9.9995E+01  1.7132E+01  4.2010E+01
 RELATIVEINF(%)  6.9994E+01  2.0780E-02  2.1437E-04  9.2236E+00  2.2429E+00
 EPSSHRINKSD(%)  3.2932E+01
 EPSSHRINKVR(%)  5.5019E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1385.7463822482173     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -650.59555568447911     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.03
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
 





 #OBJV:********************************************    -1385.746       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  4.04E-01  1.35E+00  3.20E-01  8.38E-01  1.20E+01  1.00E-02  9.32E-01  7.99E-01  2.73E+00
 


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
+        1.38E+03
 
 TH 2
+        0.00E+00  2.92E+03
 
 TH 3
+       -6.51E+00  0.00E+00  3.27E+03
 
 TH 4
+       -5.34E+01  0.00E+00 -3.42E+02  5.82E+02
 
 TH 5
+        1.55E+02  0.00E+00 -5.19E+03 -2.31E+02  9.93E+03
 
 TH 6
+       -4.08E-01  0.00E+00  1.91E+01 -1.46E+01 -6.51E+00  2.47E+02
 
 TH 7
+       -5.41E-03  0.00E+00  9.48E-03 -5.53E-02  1.03E-01  5.03E-02  1.12E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.73E+01  0.00E+00  5.74E+01 -1.26E+01  8.90E+00  5.70E+00 -1.62E-02  0.00E+00  1.57E+02
 
 TH10
+        3.79E-01  0.00E+00 -4.98E+01  7.43E+00  1.92E+01  9.22E+00 -6.92E-02  0.00E+00 -4.82E+00  1.06E+02
 
 TH11
+       -1.86E+01  0.00E+00 -1.19E+01 -8.13E+00  7.40E+00  3.92E+00 -1.70E-02  0.00E+00  6.18E+00  2.03E+01  3.84E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.461
Stop Time:
Sat Sep 18 09:59:37 CDT 2021
