Sat Sep 25 11:24:21 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.64523252141        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7984E+02 -8.3402E+01 -5.6135E+01 -6.1536E+01  5.9327E+01 -1.8151E+01 -5.3246E+00  1.5741E+01 -7.9996E+00  1.3444E+01
            -4.0322E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1626.03669899267        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9393E-01  1.0286E+00  1.2052E+00  9.6896E-01  1.0684E+00  1.0607E+00  1.0167E+00  8.6455E-01  1.0790E+00  9.2862E-01
             1.0273E+00
 PARAMETER:  9.3911E-02  1.2822E-01  2.8664E-01  6.8471E-02  1.6621E-01  1.5890E-01  1.1653E-01 -4.5551E-02  1.7601E-01  2.5948E-02
             1.2697E-01
 GRADIENT:   1.5294E+02 -6.7753E+01 -2.6539E+00 -8.9812E+01  3.2147E+01  9.5986E+00  1.9680E+00  3.0271E+00 -8.8746E-01 -1.2645E+01
             1.8608E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1629.76072525070        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.9077E-01  8.8621E-01  1.1847E+00  1.1060E+00  9.8945E-01  1.0862E+00  7.3830E-01  4.0796E-01  1.0536E+00  1.0744E+00
             1.0088E+00
 PARAMETER:  9.0724E-02 -2.0800E-02  2.6947E-01  2.0074E-01  8.9391E-02  1.8266E-01 -2.0341E-01 -7.9658E-01  1.5217E-01  1.7178E-01
             1.0875E-01
 GRADIENT:   1.4365E+02 -1.8030E+01 -9.0942E-01 -2.0866E+01  5.5174E+00  2.0368E+01 -1.3902E+00 -1.2698E+00 -9.1382E-01  4.8474E+00
            -6.7721E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.34986391587        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.2169E-01  9.0951E-01  1.0938E+00  1.1052E+00  9.6002E-01  1.0037E+00  1.1310E+00  4.3681E-01  9.6357E-01  9.8008E-01
             1.0190E+00
 PARAMETER:  1.8455E-02  5.1481E-03  1.8962E-01  1.9999E-01  5.9203E-02  1.0373E-01  2.2312E-01 -7.2826E-01  6.2885E-02  7.9884E-02
             1.1885E-01
 GRADIENT:   7.1280E-01 -6.5315E+00 -4.0687E+00  1.6753E+00  8.7784E+00 -6.8039E+00 -1.0517E-02  1.2143E-01 -2.9721E+00  4.1992E+00
             1.1469E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.75739609400        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.2334E-01  1.0371E+00  9.6179E-01  1.0223E+00  9.4249E-01  1.0191E+00  1.0609E+00  2.8548E-01  1.0205E+00  9.1285E-01
             1.0153E+00
 PARAMETER:  2.0245E-02  1.3646E-01  6.1044E-02  1.2208E-01  4.0771E-02  1.1889E-01  1.5914E-01 -1.1536E+00  1.2026E-01  8.8186E-03
             1.1515E-01
 GRADIENT:   2.9634E+00  2.2067E+00  1.7456E+00  3.1714E-01 -3.9582E+00 -7.0485E-01  1.1071E-02  1.6340E-01 -3.3763E-01  8.6562E-02
            -8.9858E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1635.10992979071        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  9.3851E-01  1.0546E+00  9.5254E-01  1.0164E+00  9.4765E-01  1.0324E+00  1.0581E+00  2.0448E-01  1.0298E+00  9.1621E-01
             1.0179E+00
 PARAMETER:  3.6542E-02  1.5317E-01  5.1381E-02  1.1624E-01  4.6229E-02  1.3185E-01  1.5650E-01 -1.4873E+00  1.2940E-01  1.2493E-02
             1.1777E-01
 GRADIENT:   5.8104E-01  1.2193E+00  7.5690E-01  9.7764E-01 -1.2064E+00 -4.1982E-02  1.4038E-01  3.1509E-02  4.1167E-02 -1.6644E-01
            -3.9882E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1635.13329788085        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  9.3864E-01  1.1127E+00  9.2038E-01  9.7739E-01  9.5986E-01  1.0330E+00  1.0181E+00  5.6069E-02  1.0581E+00  9.1852E-01
             1.0187E+00
 PARAMETER:  3.6675E-02  2.0678E-01  1.7029E-02  7.7127E-02  5.9037E-02  1.3250E-01  1.1799E-01 -2.7812E+00  1.5651E-01  1.5007E-02
             1.1849E-01
 GRADIENT:   8.0495E-02 -7.3207E-02 -5.3288E-02 -1.5054E-01 -9.2630E-03  3.1800E-02 -1.8775E-02  3.3459E-03 -8.7881E-03  3.0216E-02
            -4.9292E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1635.13494072900        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  9.3858E-01  1.1048E+00  9.2112E-01  9.8239E-01  9.5646E-01  1.0330E+00  1.0250E+00  1.4283E-02  1.0534E+00  9.1655E-01
             1.0186E+00
 PARAMETER:  3.6615E-02  1.9963E-01  1.7831E-02  8.2234E-02  5.5480E-02  1.3246E-01  1.2474E-01 -4.1487E+00  1.5200E-01  1.2865E-02
             1.1843E-01
 GRADIENT:   7.5323E-03  6.0199E-03  6.1882E-02 -1.0478E-01 -1.1371E-01  2.9445E-02  1.9087E-02  1.9430E-04 -1.0062E-02 -1.8901E-02
            -6.2513E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1635.13501389250        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  9.3858E-01  1.1058E+00  9.2116E-01  9.8178E-01  9.5706E-01  1.0329E+00  1.0238E+00  1.0000E-02  1.0541E+00  9.1710E-01
             1.0188E+00
 PARAMETER:  3.6615E-02  2.0061E-01  1.7881E-02  8.1612E-02  5.6109E-02  1.3239E-01  1.2354E-01 -4.5766E+00  1.5272E-01  1.3460E-02
             1.1858E-01
 GRADIENT:  -7.2359E-05 -3.4635E-03 -4.1168E-03  3.3734E-03  6.9992E-03  2.1191E-03 -1.5708E-03  0.0000E+00  1.3246E-03  8.1135E-04
             2.1592E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      958
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.6407E-06 -1.2897E-02 -3.3498E-04  3.6427E-03 -2.2883E-02
 SE:             2.9794E-02  1.9497E-02  1.5822E-04  2.4901E-02  2.3367E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9974E-01  5.0828E-01  3.4244E-02  8.8369E-01  3.2743E-01

 ETASHRINKSD(%)  1.8576E-01  3.4684E+01  9.9470E+01  1.6579E+01  2.1718E+01
 ETASHRINKVR(%)  3.7118E-01  5.7338E+01  9.9997E+01  3.0409E+01  3.8719E+01
 EBVSHRINKSD(%)  4.4757E-01  3.4198E+01  9.9511E+01  1.6925E+01  2.0358E+01
 EBVSHRINKVR(%)  8.9314E-01  5.6701E+01  9.9998E+01  3.0985E+01  3.6571E+01
 RELATIVEINF(%)  9.8810E+01  1.7864E+00  3.0399E-04  3.6833E+00  6.8427E+00
 EPSSHRINKSD(%)  4.2900E+01
 EPSSHRINKVR(%)  6.7396E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.1350138924984     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.98418732876019     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.135       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  1.11E+00  9.21E-01  9.82E-01  9.57E-01  1.03E+00  1.02E+00  1.00E-02  1.05E+00  9.17E-01  1.02E+00
 


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
+        1.18E+03
 
 TH 2
+       -5.58E+00  3.96E+02
 
 TH 3
+        1.16E+01  1.69E+02  3.34E+02
 
 TH 4
+       -7.88E+00  3.37E+02 -1.43E+02  7.23E+02
 
 TH 5
+       -4.43E+00 -3.16E+02 -4.52E+02  1.84E+02  9.11E+02
 
 TH 6
+       -2.87E+00 -1.94E+00  2.03E+00 -3.93E+00 -2.32E+00  1.86E+02
 
 TH 7
+        1.27E+00  1.84E+01  1.04E+01 -8.19E+00 -1.59E+01  1.25E+00  4.32E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.74E-01 -2.22E+01 -2.22E+01  2.76E+01  4.44E+00 -1.01E+00  2.24E+01  0.00E+00  9.17E+01
 
 TH10
+       -1.35E+00 -7.07E+00 -4.29E+01 -1.35E+01 -5.74E+01  1.25E+00  1.07E+01  0.00E+00  7.31E+00  1.04E+02
 
 TH11
+       -9.42E+00 -1.76E+01 -3.48E+01 -2.26E+00  1.02E+01  2.16E+00  7.00E+00  0.00E+00  8.08E+00  2.61E+01  2.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.922
Stop Time:
Sat Sep 25 11:24:39 CDT 2021
