Sat Sep 25 11:49:20 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1646.66093321596        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9121E-01 -6.2928E-01 -3.0251E+01  4.0721E+01  1.1051E+02  3.7110E+01 -6.4137E+00 -1.8253E+00 -2.9220E+01 -2.2861E+01
            -3.1987E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1656.95053415546        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0179E+00  9.2199E-01  9.5162E-01  1.0133E+00  8.7699E-01  8.9318E-01  9.8049E-01  9.9229E-01  1.1357E+00  9.9830E-01
             1.0502E+00
 PARAMETER:  1.1769E-01  1.8776E-02  5.0413E-02  1.1326E-01 -3.1262E-02 -1.2969E-02  8.0299E-02  9.2257E-02  2.2727E-01  9.8297E-02
             1.4895E-01
 GRADIENT:   4.1024E+01 -4.6139E+00  7.4045E-01 -8.4857E+00  2.1873E+00 -3.6304E+00  4.7232E+00  3.2953E+00  7.8215E+00  3.8779E+00
            -4.1001E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1657.85659814580        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0191E+00  8.3897E-01  8.6789E-01  1.0758E+00  7.9902E-01  9.0178E-01  8.5360E-01  8.0089E-01  1.1009E+00  9.3823E-01
             1.0829E+00
 PARAMETER:  1.1888E-01 -7.5585E-02 -4.1690E-02  1.7308E-01 -1.2437E-01 -3.3841E-03 -5.8293E-02 -1.2204E-01  1.9609E-01  3.6238E-02
             1.7966E-01
 GRADIENT:   3.8865E+01  8.7462E+00  3.9524E-01  2.2218E+01 -5.1951E+00 -1.8851E-01 -8.1174E-01  1.0917E+00  5.3374E+00 -7.1431E-01
             7.4480E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1658.60801652932        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0071E+00  9.8649E-01  6.2358E-01  9.5680E-01  7.3378E-01  9.1057E-01  1.1092E+00  4.7169E-01  1.0706E+00  8.1697E-01
             1.0534E+00
 PARAMETER:  1.0703E-01  8.6395E-02 -3.7228E-01  5.5839E-02 -2.0954E-01  6.3190E-03  2.0362E-01 -6.5143E-01  1.6819E-01 -1.0215E-01
             1.5200E-01
 GRADIENT:  -2.1255E+00  4.5250E+00 -5.1305E+00  6.9736E+00 -1.0471E+00  2.0088E+00  1.4833E+00  2.3410E+00  2.8063E+00  5.7038E+00
            -1.5225E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1659.10786407122        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0075E+00  8.7461E-01  5.7225E-01  1.0036E+00  6.5874E-01  9.0514E-01  1.2284E+00  2.1805E-01  9.9415E-01  7.1732E-01
             1.0609E+00
 PARAMETER:  1.0748E-01 -3.3972E-02 -4.5818E-01  1.0357E-01 -3.1743E-01  3.3371E-04  3.0569E-01 -1.4230E+00  9.4132E-02 -2.3223E-01
             1.5915E-01
 GRADIENT:  -3.9420E+00 -1.4404E+00 -7.8528E-01 -9.7756E-01  5.6669E-01 -4.4136E-01 -8.6948E-01  4.3522E-01 -1.6768E+00 -4.4953E-01
            -1.4484E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1659.47554733113        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  1.0095E+00  8.9527E-01  5.7752E-01  9.9410E-01  6.7038E-01  9.0535E-01  1.2159E+00  8.6933E-02  1.0125E+00  7.4469E-01
             1.0620E+00
 PARAMETER:  1.0943E-01 -1.0632E-02 -4.4901E-01  9.4086E-02 -2.9991E-01  5.6819E-04  2.9547E-01 -2.3426E+00  1.1244E-01 -1.9479E-01
             1.6016E-01
 GRADIENT:   2.2830E+00  4.7143E-01 -3.5777E-01 -7.6868E-01 -1.7011E+00 -1.6008E-01  6.5791E-01  7.4230E-02  5.8154E-01  1.2396E+00
             8.0757E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.07648201430        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0242E+00  8.7674E-01  6.3435E-01  1.0185E+00  7.0001E-01  9.1224E-01  1.2485E+00  4.4168E-02  1.0096E+00  8.0063E-01
             1.0609E+00
 PARAMETER:  1.2388E-01 -3.1549E-02 -3.5515E-01  1.1831E-01 -2.5666E-01  8.1485E-03  3.2195E-01 -3.0197E+00  1.0953E-01 -1.2235E-01
             1.5914E-01
 GRADIENT:  -7.8785E-01 -1.3390E+00  6.9637E-02 -2.2480E+00  9.1924E-01 -1.1681E-02 -1.0831E-01  1.0113E-02 -3.8945E-01 -1.3256E-01
            -3.3925E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.08449238520        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  1.0245E+00  8.7668E-01  6.3200E-01  1.0197E+00  6.9798E-01  9.1232E-01  1.2489E+00  1.0000E-02  1.0104E+00  7.9888E-01
             1.0617E+00
 PARAMETER:  1.2425E-01 -3.1610E-02 -3.5887E-01  1.1950E-01 -2.5956E-01  8.2354E-03  3.2225E-01 -4.5682E+00  1.1039E-01 -1.2455E-01
             1.5991E-01
 GRADIENT:   6.5267E-02 -1.5085E-02  1.1929E-02 -6.7012E-02  1.1243E-02  7.4638E-03 -1.5742E-02  0.0000E+00 -8.3568E-03 -1.5228E-02
            -8.8779E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1660.08451827395        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.0245E+00  8.7528E-01  6.3224E-01  1.0206E+00  6.9757E-01  9.1229E-01  1.2507E+00  1.0000E-02  1.0098E+00  7.9890E-01
             1.0618E+00
 PARAMETER:  1.2422E-01 -3.3210E-02 -3.5848E-01  1.2036E-01 -2.6016E-01  8.2019E-03  3.2371E-01 -4.5569E+00  1.0974E-01 -1.2452E-01
             1.5993E-01
 GRADIENT:   3.7055E-03  8.6214E-03  1.6369E-03  1.5688E-02 -5.9820E-03 -4.4340E-04 -8.0375E-04  0.0000E+00  2.9222E-03 -1.0061E-03
             2.5888E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      820
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5482E-05 -1.4872E-03 -4.5261E-04 -1.5382E-03 -9.4937E-03
 SE:             2.9804E-02  2.0606E-02  1.9471E-04  2.5963E-02  2.2770E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9932E-01  9.4246E-01  2.0096E-02  9.5275E-01  6.7672E-01

 ETASHRINKSD(%)  1.5410E-01  3.0969E+01  9.9348E+01  1.3021E+01  2.3718E+01
 ETASHRINKVR(%)  3.0797E-01  5.2347E+01  9.9996E+01  2.4347E+01  4.1811E+01
 EBVSHRINKSD(%)  5.6885E-01  3.0773E+01  9.9402E+01  1.3033E+01  2.3357E+01
 EBVSHRINKVR(%)  1.1345E+00  5.2076E+01  9.9996E+01  2.4367E+01  4.1259E+01
 RELATIVEINF(%)  9.8622E+01  3.6517E+00  4.2964E-04  8.6226E+00  4.5492E+00
 EPSSHRINKSD(%)  4.4465E+01
 EPSSHRINKVR(%)  6.9159E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.0845182739542     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -924.93369171021607     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.085       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  8.75E-01  6.32E-01  1.02E+00  6.98E-01  9.12E-01  1.25E+00  1.00E-02  1.01E+00  7.99E-01  1.06E+00
 


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
+        1.26E+03
 
 TH 2
+       -8.31E+00  4.70E+02
 
 TH 3
+        2.01E+01  3.61E+02  1.08E+03
 
 TH 4
+       -1.35E+01  2.84E+02 -3.46E+02  7.98E+02
 
 TH 5
+       -2.28E+00 -6.23E+02 -1.31E+03  3.66E+02  2.07E+03
 
 TH 6
+       -4.43E-01 -2.61E+00  4.49E+00 -2.24E+00 -1.00E+00  2.35E+02
 
 TH 7
+        1.10E+00  3.01E+01 -1.15E+01 -4.92E+00 -1.28E+01  4.56E-01  3.26E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.90E+00 -2.47E+01 -2.64E+01  2.59E+01 -2.84E+00  1.65E-01  1.49E+01  0.00E+00  1.14E+02
 
 TH10
+       -1.51E+00 -7.89E+00 -1.09E+02 -2.27E+01 -3.55E+01  4.35E-02  1.67E+01  0.00E+00  9.02E+00  1.18E+02
 
 TH11
+       -8.91E+00 -1.04E+01 -3.57E+01 -6.80E+00  1.49E+01  1.98E+00  5.57E+00  0.00E+00  6.76E+00  2.07E+01  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.502
Stop Time:
Sat Sep 25 11:49:35 CDT 2021
