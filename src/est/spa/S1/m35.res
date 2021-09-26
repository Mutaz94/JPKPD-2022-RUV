Sat Sep 25 09:50:54 CDT 2021
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
$DATA ../../../../data/spa/S1/dat35.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1693.19242683615        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.5502E+01 -3.9320E+01 -5.4499E+01  7.9594E+00  9.1804E+01  1.6285E+01 -1.2471E+00  1.1053E+01 -2.9467E+00 -1.3524E+01
            -8.2227E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1701.55625587718        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      176
 NPARAMETR:  1.0597E+00  1.0167E+00  1.0870E+00  1.0023E+00  9.7421E-01  9.2785E-01  9.9311E-01  9.2988E-01  1.0267E+00  1.0330E+00
             1.0383E+00
 PARAMETER:  1.5797E-01  1.1658E-01  1.8343E-01  1.0227E-01  7.3875E-02  2.5111E-02  9.3085E-02  2.7299E-02  1.2631E-01  1.3250E-01
             1.3756E-01
 GRADIENT:   7.6841E+00  1.6100E+00 -2.1528E+00 -3.1431E+00 -5.4705E+00 -1.1419E+01  2.4572E+00  6.3154E+00  3.2960E-01 -4.8614E+00
             6.2493E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.15873808108        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.0556E+00  9.4576E-01  1.2172E+00  1.0559E+00  1.0126E+00  9.4933E-01  8.8650E-01  7.2258E-01  1.0282E+00  1.1781E+00
             1.0201E+00
 PARAMETER:  1.5415E-01  4.4238E-02  2.9656E-01  1.5438E-01  1.1256E-01  4.7999E-02 -2.0477E-02 -2.2492E-01  1.2781E-01  2.6389E-01
             1.1994E-01
 GRADIENT:   1.0758E+00  5.7848E+00  3.5175E+00  5.8862E+00 -1.2052E+00 -1.5290E+00  4.4050E-02 -1.1217E+00 -2.9920E-02  7.1411E-01
             4.8624E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1703.38583595091        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0544E+00  8.0827E-01  1.1807E+00  1.1362E+00  9.4230E-01  9.5317E-01  1.0306E+00  6.9614E-01  9.5274E-01  1.1055E+00
             1.0202E+00
 PARAMETER:  1.5294E-01 -1.1286E-01  2.6612E-01  2.2765E-01  4.0573E-02  5.2034E-02  1.3010E-01 -2.6221E-01  5.1586E-02  2.0025E-01
             1.2000E-01
 GRADIENT:  -2.8039E-01  1.1357E+00  7.4205E-01 -4.0574E-01 -1.8159E+00  1.7977E-01 -2.5822E-01  7.3761E-02 -2.9064E-01 -9.8324E-01
             6.4773E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1703.57921218361        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0510E+00  5.6961E-01  1.3429E+00  1.2953E+00  9.2365E-01  9.4824E-01  1.1722E+00  8.0526E-01  8.6974E-01  1.1393E+00
             1.0180E+00
 PARAMETER:  1.4972E-01 -4.6281E-01  3.9487E-01  3.5871E-01  2.0581E-02  4.6854E-02  2.5890E-01 -1.1659E-01 -3.9561E-02  2.3038E-01
             1.1783E-01
 GRADIENT:  -5.8420E-01  3.9963E+00  2.1622E+00  1.0383E+01 -2.6990E+00 -5.9980E-01 -8.7199E-02 -4.5334E-01 -7.0363E-01  1.0800E-01
            -1.7794E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.69685578564        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0492E+00  4.3331E-01  1.3654E+00  1.3760E+00  8.9012E-01  9.4919E-01  1.3590E+00  8.3418E-01  8.2845E-01  1.1269E+00
             1.0168E+00
 PARAMETER:  1.4807E-01 -7.3630E-01  4.1145E-01  4.1921E-01 -1.6401E-02  4.7851E-02  4.0676E-01 -8.1303E-02 -8.8200E-02  2.1944E-01
             1.1666E-01
 GRADIENT:  -1.6862E-01  1.6866E+00  1.2399E+00  5.9272E+00 -3.2188E+00  4.0545E-01 -1.2239E-01 -1.5976E-01  6.7299E-02  6.0088E-01
            -4.7691E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.75437429638        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0474E+00  3.1285E-01  1.4137E+00  1.4501E+00  8.7304E-01  9.4746E-01  1.7525E+00  8.9312E-01  7.8994E-01  1.1135E+00
             1.0181E+00
 PARAMETER:  1.4626E-01 -1.0620E+00  4.4622E-01  4.7164E-01 -3.5779E-02  4.6034E-02  6.6102E-01 -1.3030E-02 -1.3580E-01  2.0752E-01
             1.1798E-01
 GRADIENT:   7.5941E-03  8.5024E-01  1.0072E+00  4.3964E+00 -1.3928E+00  3.2232E-01  5.6457E-02 -1.4998E-01  1.5135E-01 -3.3671E-01
            -9.9696E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.77867225877        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  1.0456E+00  2.1827E-01  1.4424E+00  1.5078E+00  8.5689E-01  9.4516E-01  2.2371E+00  9.3199E-01  7.6427E-01  1.1118E+00
             1.0183E+00
 PARAMETER:  1.4461E-01 -1.4220E+00  4.6628E-01  5.1066E-01 -5.4449E-02  4.3595E-02  9.0518E-01  2.9568E-02 -1.6884E-01  2.0601E-01
             1.1810E-01
 GRADIENT:  -2.2525E-01  6.3860E-01  4.9888E-01  4.3561E+00 -2.1085E+00 -1.2361E-01  1.0697E-01  4.5741E-02 -5.3965E-03  5.5743E-01
            -2.7002E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.82426885884        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0440E+00  1.2204E-01  1.5490E+00  1.5707E+00  8.6743E-01  9.4316E-01  2.6145E+00  1.0317E+00  7.4434E-01  1.1299E+00
             1.0179E+00
 PARAMETER:  1.4303E-01 -2.0034E+00  5.3761E-01  5.5155E-01 -4.2215E-02  4.1478E-02  1.0611E+00  1.3117E-01 -1.9526E-01  2.2213E-01
             1.1769E-01
 GRADIENT:   1.3224E-03  2.0944E-01 -2.3030E-01  2.8733E+00  6.2242E-01 -3.4216E-01 -1.1318E-02 -9.4763E-02 -4.5232E-01  3.2692E-01
            -2.9805E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.86244274335        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1590
 NPARAMETR:  1.0426E+00  5.5842E-02  1.5829E+00  1.6121E+00  8.5790E-01  9.4307E-01  2.8114E+00  1.0761E+00  7.3340E-01  1.1221E+00
             1.0195E+00
 PARAMETER:  1.4174E-01 -2.7852E+00  5.5924E-01  5.7752E-01 -5.3268E-02  4.1389E-02  1.1337E+00  1.7330E-01 -2.1006E-01  2.1523E-01
             1.1929E-01
 GRADIENT:  -4.2584E-01  1.5894E-01  9.4014E-01  3.9745E+00 -2.5571E+00 -8.4336E-03 -6.7803E-03  5.0050E-02  4.3483E-01  1.9924E-02
             1.8550E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1703.89571762338        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1766
 NPARAMETR:  1.0421E+00  2.0218E-02  1.6118E+00  1.6329E+00  8.5981E-01  9.4217E-01  2.7087E+00  1.1012E+00  7.2495E-01  1.1273E+00
             1.0189E+00
 PARAMETER:  1.4119E-01 -3.8012E+00  5.7738E-01  5.9033E-01 -5.1047E-02  4.0427E-02  1.0965E+00  1.9644E-01 -2.2166E-01  2.1983E-01
             1.1872E-01
 GRADIENT:  -1.8228E-01  2.6929E-03 -1.4226E-01 -1.5276E+00  2.7760E-01 -1.7831E-01  9.7768E-04  3.2645E-02  4.9673E-01  1.8956E-01
             1.0468E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1703.90523909077        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1941
 NPARAMETR:  1.0420E+00  1.0000E-02  1.6191E+00  1.6400E+00  8.5906E-01  9.4258E-01  2.5784E+00  1.1092E+00  7.2108E-01  1.1251E+00
             1.0188E+00
 PARAMETER:  1.4112E-01 -4.5162E+00  5.8188E-01  5.9470E-01 -5.1917E-02  4.0864E-02  1.0472E+00  2.0367E-01 -2.2700E-01  2.1787E-01
             1.1865E-01
 GRADIENT:   5.3779E-02  0.0000E+00 -3.9057E-02 -7.4279E-03  1.0062E-01  5.0493E-02  2.3969E-04 -1.4174E-03 -3.4543E-02 -2.5129E-02
            -1.1655E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1703.90524982352        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2100
 NPARAMETR:  1.0419E+00  1.0000E-02  1.6189E+00  1.6400E+00  8.5893E-01  9.4243E-01  2.5778E+00  1.1089E+00  7.2115E-01  1.1252E+00
             1.0188E+00
 PARAMETER:  1.4109E-01 -4.5163E+00  5.8175E-01  5.9469E-01 -5.2068E-02  4.0705E-02  1.0469E+00  2.0358E-01 -2.2691E-01  2.1792E-01
             1.1867E-01
 GRADIENT:  -6.1211E-02  0.0000E+00 -3.1149E-03  2.4612E-01 -9.5207E-03 -1.7513E-02  1.6182E-03  4.4302E-03 -7.5297E-03 -1.0181E-03
             5.5834E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2100
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2056E-04 -8.7080E-04 -2.4747E-02 -6.8601E-03 -3.2990E-02
 SE:             2.9852E-02  4.5806E-04  1.4774E-02  2.9187E-02  2.2627E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9678E-01  5.7295E-02  9.3917E-02  8.1418E-01  1.4484E-01

 ETASHRINKSD(%)  1.0000E-10  9.8465E+01  5.0507E+01  2.2198E+00  2.4197E+01
 ETASHRINKVR(%)  1.0000E-10  9.9976E+01  7.5504E+01  4.3903E+00  4.2539E+01
 EBVSHRINKSD(%)  4.6055E-01  9.8570E+01  5.3978E+01  2.5098E+00  2.0244E+01
 EBVSHRINKVR(%)  9.1899E-01  9.9980E+01  7.8820E+01  4.9565E+00  3.6390E+01
 RELATIVEINF(%)  9.6486E+01  9.8024E-04  4.3664E+00  5.3474E+00  9.0608E+00
 EPSSHRINKSD(%)  4.4134E+01
 EPSSHRINKVR(%)  6.8790E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.9052498235224     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.75442325978418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.905       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.62E+00  1.64E+00  8.59E-01  9.42E-01  2.58E+00  1.11E+00  7.21E-01  1.13E+00  1.02E+00
 


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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.96E+00  0.00E+00  7.86E+01
 
 TH 4
+       -1.17E+01  0.00E+00 -3.07E+01  7.55E+02
 
 TH 5
+        2.48E+01  0.00E+00 -1.94E+02 -5.65E+01  7.70E+02
 
 TH 6
+        2.28E+01  0.00E+00  3.57E-02  1.14E+00 -2.17E+00  2.16E+02
 
 TH 7
+       -6.47E-01  0.00E+00 -9.05E-02  2.06E-02  1.91E+00 -5.57E-01 -6.78E-03
 
 TH 8
+       -6.25E+00  0.00E+00 -2.05E+01 -1.66E+00 -1.88E+01 -2.98E+00  1.62E-01  1.94E+01
 
 TH 9
+        4.03E+00  0.00E+00  6.14E+00 -2.71E+00 -9.51E+00 -2.57E+00  2.21E-02  1.03E+00  3.48E+02
 
 TH10
+       -5.49E+00  0.00E+00 -3.96E+00 -2.32E+00 -7.71E+01  4.25E+00 -3.63E-01  1.70E+01 -1.49E+00  5.85E+01
 
 TH11
+        1.06E+01  0.00E+00 -1.09E+01 -9.69E+00 -3.38E+01 -3.88E+01 -3.94E-01  5.67E+00  1.08E+01  1.45E+01  2.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.130
Stop Time:
Sat Sep 25 09:51:28 CDT 2021
