Sat Sep 18 08:40:59 CDT 2021
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
$DATA ../../../../data/spa/B/dat70.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1726.57930250076        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.6190E+01 -2.1083E+01 -9.9462E+00 -3.1902E+01 -2.3949E+01  5.7114E+00 -2.9607E+00  1.1273E+01  1.2662E+00  2.2599E+01
             4.5724E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1733.01371980271        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0316E+00  1.0196E+00  1.0560E+00  1.0044E+00  1.0492E+00  9.7811E-01  1.0338E+00  9.4310E-01  9.9441E-01  9.1470E-01
             8.8181E-01
 PARAMETER:  1.3112E-01  1.1944E-01  1.5452E-01  1.0443E-01  1.4805E-01  7.7864E-02  1.3321E-01  4.1412E-02  9.4398E-02  1.0844E-02
            -2.5779E-02
 GRADIENT:   7.5636E+01 -4.5100E+00  1.0823E+00 -5.5100E+00  8.8653E+00 -4.5658E-01 -4.6026E+00  2.8516E+00 -2.1269E-01 -6.8867E-01
            -6.7372E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1736.03341880976        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0251E+00  8.7042E-01  7.8719E-01  1.0885E+00  8.3711E-01  1.0055E+00  1.4451E+00  5.0950E-01  8.7528E-01  6.9545E-01
             8.7710E-01
 PARAMETER:  1.2477E-01 -3.8782E-02 -1.3929E-01  1.8481E-01 -7.7795E-02  1.0545E-01  4.6820E-01 -5.7433E-01 -3.3210E-02 -2.6320E-01
            -3.1134E-02
 GRADIENT:   5.3069E+01  7.6171E+00 -1.9578E+01  4.3473E+01  2.7985E+01  9.6437E+00  9.4815E+00  2.8377E+00  5.2971E+00 -3.6030E-01
            -7.3573E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1736.35592722825        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0167E+00  9.0152E-01  7.4772E-01  1.0563E+00  8.2518E-01  9.9640E-01  1.3523E+00  4.2156E-01  8.8073E-01  6.9417E-01
             8.8546E-01
 PARAMETER:  1.1654E-01 -3.6734E-03 -1.9073E-01  1.5473E-01 -9.2156E-02  9.6394E-02  4.0178E-01 -7.6378E-01 -2.7003E-02 -2.6503E-01
            -2.1649E-02
 GRADIENT:   2.6995E+01  1.0371E+00 -1.2300E+01  2.0290E+01  1.5189E+01  5.1349E+00  3.8743E+00  2.0423E+00  2.9879E+00  1.8152E+00
            -2.9987E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1736.36867904824        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.0137E+00  9.1174E-01  7.1650E-01  1.0435E+00  8.0856E-01  9.9284E-01  1.3299E+00  3.5783E-01  8.8040E-01  6.7114E-01
             8.8863E-01
 PARAMETER:  1.1357E-01  7.5942E-03 -2.3338E-01  1.4258E-01 -1.1250E-01  9.2816E-02  3.8512E-01 -9.2770E-01 -2.7374E-02 -2.9877E-01
            -1.8069E-02
 GRADIENT:   1.6914E+01  3.3297E-01 -8.4594E+00  1.2734E+01  9.8003E+00  3.2994E+00  2.4698E+00  1.4965E+00  1.9141E+00  1.5188E+00
            -1.7066E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1737.38141857377        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  1.0351E+00  7.9303E-01  7.1974E-01  1.1132E+00  7.6016E-01  9.9806E-01  1.5009E+00  2.1167E-01  8.3073E-01  6.6353E-01
             8.9571E-01
 PARAMETER:  1.3452E-01 -1.3190E-01 -2.2887E-01  2.0722E-01 -1.7422E-01  9.8055E-02  5.0603E-01 -1.4527E+00 -8.5449E-02 -3.1019E-01
            -1.0133E-02
 GRADIENT:   2.3013E+00  1.9662E-02 -1.7574E+00  1.7405E+00  7.3528E-01  3.6685E-01  7.7948E-02  2.2237E-01  1.2268E-01  5.7574E-01
            -9.8022E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1737.45085327407        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  1.0348E+00  8.1272E-01  7.0438E-01  1.0986E+00  7.5826E-01  9.9741E-01  1.4653E+00  9.0113E-02  8.3561E-01  6.6358E-01
             8.9764E-01
 PARAMETER:  1.3417E-01 -1.0737E-01 -2.5043E-01  1.9399E-01 -1.7672E-01  9.7407E-02  4.8204E-01 -2.3067E+00 -7.9588E-02 -3.1011E-01
            -7.9893E-03
 GRADIENT:   8.1609E-01 -2.4677E-01  8.5788E-01 -1.8876E+00 -2.7192E+00 -1.5144E-03 -5.5021E-01  3.4549E-02 -3.9258E-01  7.5302E-01
             7.7217E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1737.47217461993        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.0345E+00  8.2484E-01  7.0212E-01  1.0920E+00  7.6256E-01  9.9749E-01  1.4534E+00  1.9970E-02  8.4168E-01  6.6157E-01
             8.9601E-01
 PARAMETER:  1.3389E-01 -9.2563E-02 -2.5365E-01  1.8806E-01 -1.7107E-01  9.7483E-02  4.7392E-01 -3.8135E+00 -7.2354E-02 -3.1313E-01
            -9.8017E-03
 GRADIENT:   1.3433E-02 -9.8064E-02 -2.0273E-01  2.2319E-01  3.8576E-01 -2.3650E-02 -7.2930E-02  1.4489E-03 -6.0612E-03  5.7314E-03
            -4.5345E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1737.47283082467        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  1.0345E+00  8.2367E-01  7.0132E-01  1.0925E+00  7.6146E-01  9.9755E-01  1.4558E+00  1.0000E-02  8.4120E-01  6.6036E-01
             8.9614E-01
 PARAMETER:  1.3389E-01 -9.3983E-02 -2.5479E-01  1.8846E-01 -1.7251E-01  9.7546E-02  4.7558E-01 -4.5527E+00 -7.2926E-02 -3.1497E-01
            -9.6551E-03
 GRADIENT:  -1.6692E-03 -9.7341E-03 -8.1774E-03 -1.0896E-02  1.0761E-02  1.0494E-03  1.9751E-03  0.0000E+00  1.5747E-03  2.0368E-03
             2.3387E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      965
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1772E-04  9.7576E-03 -5.7788E-04 -9.3557E-03 -1.7584E-03
 SE:             2.9877E-02  2.2541E-02  2.3766E-04  2.5296E-02  2.1258E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9686E-01  6.6511E-01  1.5034E-02  7.1150E-01  9.3408E-01

 ETASHRINKSD(%)  1.0000E-10  2.4483E+01  9.9204E+01  1.5255E+01  2.8784E+01
 ETASHRINKVR(%)  1.0000E-10  4.2972E+01  9.9994E+01  2.8182E+01  4.9283E+01
 EBVSHRINKSD(%)  3.4845E-01  2.3864E+01  9.9266E+01  1.5333E+01  2.8508E+01
 EBVSHRINKVR(%)  6.9569E-01  4.2033E+01  9.9995E+01  2.8315E+01  4.8889E+01
 RELATIVEINF(%)  9.9016E+01  6.0723E+00  5.6263E-04  9.7028E+00  3.6664E+00
 EPSSHRINKSD(%)  4.4275E+01
 EPSSHRINKVR(%)  6.8947E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1737.4728308246686     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1002.3220042609304     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1737.473       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.24E-01  7.01E-01  1.09E+00  7.61E-01  9.98E-01  1.46E+00  1.00E-02  8.41E-01  6.60E-01  8.96E-01
 


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
+        1.04E+03
 
 TH 2
+       -8.60E+00  4.68E+02
 
 TH 3
+        2.10E+01  3.68E+02  1.33E+03
 
 TH 4
+       -7.94E+00  2.80E+02 -4.80E+02  9.58E+02
 
 TH 5
+       -6.99E+00 -5.86E+02 -1.54E+03  5.34E+02  2.24E+03
 
 TH 6
+       -4.79E-01 -7.15E-01  3.96E+00 -1.89E+00 -1.58E+00  1.99E+02
 
 TH 7
+        1.02E+00  3.76E+01 -2.11E+01 -1.22E+01  6.76E+00 -8.81E-02  3.71E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.37E-02 -2.52E+01 -4.02E+01  3.00E+01  1.65E-01 -2.02E+00  1.22E+01  0.00E+00  1.48E+02
 
 TH10
+       -1.03E+00 -1.44E+01 -1.22E+02 -3.71E+01 -4.54E+01 -3.94E-01  1.36E+01  0.00E+00  2.19E+01  1.33E+02
 
 TH11
+       -7.33E+00 -1.02E+01 -5.09E+01 -7.41E+00  1.87E+01  1.06E+00  4.91E+00  0.00E+00  1.47E+01  3.01E+01  2.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.918
Stop Time:
Sat Sep 18 08:41:17 CDT 2021
