Thu Sep 30 03:44:48 CDT 2021
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
$DATA ../../../../data/spa1/D/dat88.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   27392.7899750036        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.1868E+02  5.8970E+02 -1.6951E+01  5.8095E+02  2.6511E+02 -2.8420E+03 -1.2733E+03 -9.4982E+01 -1.8627E+03 -7.0283E+02
            -5.1509E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -449.697811598304        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2128E+00  1.0593E+00  8.7870E-01  1.6893E+00  1.4334E+00  2.3711E+00  1.2251E+00  9.4858E-01  1.3374E+00  9.7631E-01
             1.4168E+01
 PARAMETER:  2.9295E-01  1.5760E-01 -2.9313E-02  6.2429E-01  4.6005E-01  9.6337E-01  3.0305E-01  4.7210E-02  3.9070E-01  7.6029E-02
             2.7510E+00
 GRADIENT:  -1.2891E+01  5.2313E+01 -6.4690E+00  7.8030E+01 -5.0864E+00  6.2994E+01 -8.5447E+00  5.0948E+00 -4.4904E+01  1.0411E+00
            -8.5937E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -487.022654388314        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2416E+00  8.7732E-01  9.3263E-01  1.9339E+00  4.9598E+00  2.0774E+00  1.0279E+00  3.0094E-01  2.5603E+00  6.3533E-01
             1.3955E+01
 PARAMETER:  3.1640E-01 -3.0879E-02  3.0254E-02  7.5952E-01  1.7014E+00  8.3112E-01  1.2747E-01 -1.1008E+00  1.0401E+00 -3.5361E-01
             2.7358E+00
 GRADIENT:  -1.0164E+01  2.9682E+01 -4.9928E+00  6.1120E+01 -1.3510E+01 -3.3102E+01  3.6392E+00  3.6009E-01  2.4189E+01  2.5083E-01
             7.0112E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -535.656858218718        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0496E+00  2.4023E-01  3.6496E-01  1.5064E+00  1.7468E+01  2.0240E+00  1.7745E-01  1.0000E-02  1.0647E+00  8.9710E+00
             1.3365E+01
 PARAMETER:  1.4845E-01 -1.3261E+00 -9.0797E-01  5.0969E-01  2.9604E+00  8.0509E-01 -1.6291E+00 -6.5260E+00  1.6268E-01  2.2940E+00
             2.6927E+00
 GRADIENT:  -2.4589E+00  1.2932E+01  1.7793E+01  6.8534E+01  2.6549E+00  1.1796E+01  4.9509E-02  0.0000E+00 -3.0516E+01 -9.0891E-02
             2.3249E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -584.184749151678        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  7.6654E-01  4.6451E-02  5.6000E-02  6.1855E-01  1.7397E+02  1.8336E+00  1.0000E-02  1.0000E-02  3.4913E-01  1.3013E+01
             1.3198E+01
 PARAMETER: -1.6587E-01 -2.9694E+00 -2.7824E+00 -3.8038E-01  5.2589E+00  7.0631E-01 -4.6270E+00 -1.5228E+01 -9.5230E-01  2.6659E+00
             2.6801E+00
 GRADIENT:   1.3816E+02 -1.3187E+01 -1.1430E+02  2.1067E+02  1.0699E-01 -4.4306E+00  0.0000E+00  0.0000E+00 -2.6135E+00  3.4845E-02
            -4.0474E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -587.199562116470        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      423
 NPARAMETR:  6.9931E-01  3.5851E-02  4.1701E-02  5.0219E-01  2.6867E+02  1.8146E+00  1.0000E-02  1.0000E-02  2.9862E-01  1.2567E+01
             1.3077E+01
 PARAMETER: -2.5767E-01 -3.2284E+00 -3.0772E+00 -5.8878E-01  5.6935E+00  6.9585E-01 -5.1697E+00 -1.6516E+01 -1.1086E+00  2.6311E+00
             2.6709E+00
 GRADIENT:   1.4351E+02 -1.2424E+01 -1.7562E+02  2.0544E+02  6.3914E-02 -1.4650E+01  0.0000E+00  0.0000E+00 -1.0604E+00  1.2292E-02
            -6.8257E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -625.008005119760        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      602
 NPARAMETR:  4.0666E-01  1.0000E-02  1.9140E-02  2.4345E-01  1.5000E+03  1.6044E+00  1.0000E-02  1.0000E-02  2.4960E-01  9.1112E+00
             1.2647E+01
 PARAMETER: -7.9978E-01 -4.5957E+00 -3.8560E+00 -1.3129E+00  7.4132E+00  5.7274E-01 -8.4025E+00 -1.7316E+01 -1.2879E+00  2.3095E+00
             2.6374E+00
 GRADIENT:  -1.7328E+00  0.0000E+00  1.5943E-01  2.3910E+00  2.9730E-04 -2.3267E+00  0.0000E+00  0.0000E+00 -2.5291E-01 -1.2999E-06
            -2.5633E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -625.590136111968        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  3.8864E-01  1.0000E-02  1.6759E-02  2.1831E-01  1.7863E+03  1.6143E+00  1.0000E-02  1.0000E-02  3.4387E-01  9.1981E+00
             1.2602E+01
 PARAMETER: -8.4511E-01 -4.7466E+00 -3.9888E+00 -1.4219E+00  7.5879E+00  5.7890E-01 -8.7680E+00 -1.7868E+01 -9.6750E-01  2.3190E+00
             2.6339E+00
 GRADIENT:   7.9370E-01  0.0000E+00 -9.1632E+00  8.7700E+00  2.2987E-04  1.0871E+00  0.0000E+00  0.0000E+00  7.1287E-02 -2.4348E-07
             1.6321E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -625.815926953386        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  3.9227E-01  1.0000E-02  1.7172E-02  2.1483E-01  1.6441E+03  1.6173E+00  1.0000E-02  1.0000E-02  6.3651E-01  4.4375E+00
             1.2132E+01
 PARAMETER: -8.3581E-01 -4.7466E+00 -3.9645E+00 -1.4379E+00  7.5050E+00  5.8077E-01 -8.7680E+00 -1.7868E+01 -3.5176E-01  1.5901E+00
             2.5959E+00
 GRADIENT:  -3.9316E+00  0.0000E+00  3.5742E+00 -7.4059E+00 -1.1034E-04 -3.8433E+00  0.0000E+00  0.0000E+00  1.5419E-01  5.0366E-08
             3.0953E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -625.867950687635        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  3.9355E-01  1.0000E-02  1.7166E-02  2.1557E-01  7.7150E+01  1.6381E+00  1.0000E-02  1.0000E-02  6.4015E-01  3.6164E+00
             1.2077E+01
 PARAMETER: -8.3256E-01 -4.7466E+00 -3.9648E+00 -1.4345E+00  4.4458E+00  5.9354E-01 -8.7680E+00 -1.7868E+01 -3.4605E-01  1.3855E+00
             2.5913E+00
 GRADIENT:   5.7220E+01  0.0000E+00  7.9983E+01  2.5901E+01  6.6459E-04  9.9763E+00  0.0000E+00  0.0000E+00  3.0881E-01  8.0382E-05
             2.4817E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -625.869977880662        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:     1253             RESET HESSIAN, TYPE I
 NPARAMETR:  3.9415E-01  1.0000E-02  1.7137E-02  2.1515E-01  4.4975E+01  1.6373E+00  1.0000E-02  1.0000E-02  6.3390E-01  2.9831E+00
             1.2096E+01
 PARAMETER: -8.3103E-01 -4.7466E+00 -3.9665E+00 -1.4364E+00  3.9061E+00  5.9308E-01 -8.7680E+00 -1.7868E+01 -3.5587E-01  1.1930E+00
             2.5929E+00
 GRADIENT:   5.8292E+01  0.0000E+00  8.1370E+01  2.3654E+01  1.5142E-04  9.9998E+00  0.0000E+00  0.0000E+00  1.3839E-01  1.7450E-04
             2.5267E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -625.871209203938        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1445             RESET HESSIAN, TYPE I
 NPARAMETR:  3.9424E-01  1.0000E-02  1.7129E-02  2.1527E-01  3.9837E+01  1.6368E+00  1.0000E-02  1.0000E-02  6.3602E-01  2.6861E+00
             1.2100E+01
 PARAMETER: -8.3080E-01 -4.7466E+00 -3.9670E+00 -1.4358E+00  3.7848E+00  5.9274E-01 -8.7680E+00 -1.7868E+01 -3.5253E-01  1.0881E+00
             2.5932E+00
 GRADIENT:   5.8305E+01  0.0000E+00  7.9932E+01  2.5408E+01  1.7880E-04  9.8398E+00  0.0000E+00  0.0000E+00  2.9414E-01  1.8337E-04
             2.5509E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -625.871296238497        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1638
 NPARAMETR:  3.9435E-01  1.0000E-02  1.7126E-02  2.1527E-01  4.0823E+01  1.6366E+00  1.0000E-02  1.0000E-02  6.3494E-01  2.3244E+00
             1.2104E+01
 PARAMETER: -8.3052E-01 -4.7466E+00 -3.9671E+00 -1.4358E+00  3.8092E+00  5.9265E-01 -8.7680E+00 -1.7868E+01 -3.5423E-01  9.4345E-01
             2.5935E+00
 GRADIENT:   7.6649E-01  0.0000E+00 -4.2492E+00 -5.4485E-01 -7.9419E-05  1.3539E-01  0.0000E+00  0.0000E+00  1.6374E-02  7.6216E-05
            -3.5543E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -625.871313136674        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1835
 NPARAMETR:  3.9435E-01  1.0000E-02  1.7127E-02  2.1529E-01  4.1036E+01  1.6366E+00  1.0000E-02  1.0000E-02  6.3462E-01  1.8000E+00
             1.2105E+01
 PARAMETER: -8.3052E-01 -4.7466E+00 -3.9671E+00 -1.4358E+00  3.8144E+00  5.9259E-01 -8.7680E+00 -1.7868E+01 -3.5474E-01  6.8779E-01
             2.5936E+00
 GRADIENT:   7.5480E-01  0.0000E+00 -4.2346E+00 -5.5048E-01 -8.9739E-06  1.2195E-01  0.0000E+00  0.0000E+00  1.2884E-02  4.4337E-05
            -3.2988E-01

0ITERATION NO.:   68    OBJECTIVE VALUE:  -625.871313765519        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1949
 NPARAMETR:  3.9439E-01  1.0000E-02  1.7123E-02  2.1530E-01  4.0745E+01  1.6365E+00  1.0000E-02  1.0000E-02  6.3403E-01  1.7457E+00
             1.2106E+01
 PARAMETER: -8.3051E-01 -4.7466E+00 -3.9671E+00 -1.4358E+00  3.8143E+00  5.9259E-01 -8.7680E+00 -1.7868E+01 -3.5480E-01  6.6379E-01
             2.5936E+00
 GRADIENT:  -5.3151E-02  0.0000E+00  1.8173E-01 -7.6039E-02  3.4897E-05  9.5353E-03  0.0000E+00  0.0000E+00  1.2455E-02  3.0393E-05
            -6.2619E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1949
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6847E-03  2.3223E-06  1.8727E-04 -1.9885E-02 -1.4196E-05
 SE:             2.8970E-02  1.4837E-06  2.7619E-04  1.9948E-02  4.2606E-05
 N:                     100         100         100         100         100

 P VAL.:         9.2616E-01  1.1755E-01  4.9775E-01  3.1884E-01  7.3899E-01

 ETASHRINKSD(%)  2.9468E+00  9.9995E+01  9.9075E+01  3.3173E+01  9.9857E+01
 ETASHRINKVR(%)  5.8068E+00  1.0000E+02  9.9991E+01  5.5342E+01  1.0000E+02
 EBVSHRINKSD(%)  2.9745E+00  9.9990E+01  9.9144E+01  3.5894E+01  9.9855E+01
 EBVSHRINKVR(%)  5.8605E+00  1.0000E+02  9.9993E+01  5.8905E+01  1.0000E+02
 RELATIVEINF(%)  1.7236E+00  2.4866E-07  3.2526E-05  1.6065E-01  5.6443E-06
 EPSSHRINKSD(%)  7.5871E+00
 EPSSHRINKVR(%)  1.4599E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -625.87131376551940     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       293.06721943915329     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -625.871       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.94E-01  1.00E-02  1.71E-02  2.15E-01  4.10E+01  1.64E+00  1.00E-02  1.00E-02  6.35E-01  1.76E+00  1.21E+01
 


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
+        2.48E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.23E+04  0.00E+00  2.23E+06
 
 TH 4
+       -4.29E+01  0.00E+00 -1.98E+05  1.96E+04
 
 TH 5
+        2.75E-02  0.00E+00 -5.71E-01  3.83E-02  2.08E-06
 
 TH 6
+        4.60E+00  0.00E+00 -9.00E+00 -4.08E+01  3.59E-04  6.43E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.38E+01  0.00E+00 -2.39E+03  2.49E+02 -3.28E-04 -4.21E+00  0.00E+00  0.00E+00  2.60E+01
 
 TH10
+       -3.43E-04  0.00E+00 -4.50E-03  5.33E-04 -3.60E-07  1.03E-04  0.00E+00  0.00E+00  2.14E-03  2.37E-04
 
 TH11
+       -2.26E+01  0.00E+00  6.04E+02 -4.07E+01 -3.43E-04  1.31E+00  0.00E+00  0.00E+00  4.39E+00  1.74E-06  2.93E+00
 
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
 #CPUT: Total CPU Time in Seconds,       49.642
Stop Time:
Thu Sep 30 03:45:39 CDT 2021
