Wed Sep 29 11:06:42 CDT 2021
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
$DATA ../../../../data/spa/B/dat29.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1617.83896621359        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6964E+02 -2.6254E+01 -4.8392E+01  2.6446E+01  6.8264E+01  2.3763E+01 -2.0145E+01  6.9651E+00 -2.3550E+01  1.6881E+01
            -3.4979E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1628.58928821659        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8536E-01  1.0421E+00  1.0944E+00  1.0137E+00  1.0116E+00  1.0871E+00  1.0950E+00  9.7454E-01  1.1309E+00  9.1473E-01
             1.1130E+00
 PARAMETER:  8.5248E-02  1.4125E-01  1.9024E-01  1.1362E-01  1.1158E-01  1.8348E-01  1.9077E-01  7.4211E-02  2.2298E-01  1.0873E-02
             2.0702E-01
 GRADIENT:  -1.2886E+00 -1.0509E+01 -1.3335E+01  1.0415E+00  1.4525E+01  7.9709E+00 -6.7319E+00  2.0828E+00  2.4452E-01  4.9161E+00
             1.0750E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1629.31411031943        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8824E-01  9.2539E-01  1.0711E+00  1.0862E+00  9.3790E-01  1.0516E+00  1.3529E+00  9.3411E-01  1.0664E+00  8.0804E-01
             1.0806E+00
 PARAMETER:  8.8175E-02  2.2465E-02  1.6868E-01  1.8267E-01  3.5886E-02  1.5027E-01  4.0222E-01  3.1842E-02  1.6431E-01 -1.1314E-01
             1.7749E-01
 GRADIENT:   6.0344E+00 -6.1441E-02 -4.8223E+00  9.8188E-01 -7.2744E-01 -5.2036E+00  2.8159E+00  1.8971E+00  4.7653E+00 -8.1105E-01
            -1.6766E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1629.73221103960        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8368E-01  8.4724E-01  1.2345E+00  1.1490E+00  9.7195E-01  1.0653E+00  1.3505E+00  9.9765E-01  1.0297E+00  8.7625E-01
             1.0881E+00
 PARAMETER:  8.3547E-02 -6.5773E-02  3.1065E-01  2.3893E-01  7.1547E-02  1.6326E-01  4.0047E-01  9.7647E-02  1.2928E-01 -3.2100E-02
             1.8442E-01
 GRADIENT:  -7.0451E-01  5.4164E+00  1.9730E+00  5.3741E+00 -3.9061E+00  4.9888E-01 -1.6875E-01 -2.2198E-01 -4.2981E-01  1.6118E-01
            -1.8071E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1629.93815967082        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.8106E-01  5.6748E-01  1.4873E+00  1.3401E+00  9.6170E-01  1.0607E+00  1.5698E+00  1.1647E+00  9.4916E-01  8.9522E-01
             1.0910E+00
 PARAMETER:  8.0873E-02 -4.6656E-01  4.9699E-01  3.9272E-01  6.0951E-02  1.5893E-01  5.5098E-01  2.5249E-01  4.7822E-02 -1.0691E-02
             1.8707E-01
 GRADIENT:   4.6957E-01  6.7740E+00  4.3915E+00  1.3039E+01 -6.0605E+00 -1.3339E-01 -6.3588E-01 -5.3712E-01 -1.4646E+00 -8.0972E-01
             1.9752E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1630.31738967502        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.7761E-01  3.2178E-01  1.6070E+00  1.4948E+00  9.3133E-01  1.0579E+00  1.9843E+00  1.2360E+00  8.8978E-01  9.0431E-01
             1.0881E+00
 PARAMETER:  7.7354E-02 -1.0339E+00  5.7435E-01  5.0201E-01  2.8853E-02  1.5629E-01  7.8529E-01  3.1191E-01 -1.6781E-02 -5.8455E-04
             1.8440E-01
 GRADIENT:  -2.2167E-01  3.0394E+00  1.1362E+00  7.5478E+00 -3.0017E+00 -3.6288E-01 -6.0583E-01 -4.8346E-01 -7.3599E-01  1.7909E-01
            -6.3056E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1630.33215697444        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7659E-01  2.5444E-01  1.6625E+00  1.5428E+00  9.2754E-01  1.0574E+00  2.1935E+00  1.2841E+00  8.7325E-01  9.0691E-01
             1.0882E+00
 PARAMETER:  7.6314E-02 -1.2687E+00  6.0832E-01  5.3360E-01  2.4776E-02  1.5584E-01  8.8552E-01  3.5003E-01 -3.5529E-02  2.2863E-03
             1.8456E-01
 GRADIENT:  -5.1195E-01  3.3666E+00  1.6252E+00  1.2964E+01 -5.4064E+00 -3.4596E-01 -5.3232E-01 -2.4373E-01 -8.3143E-01  6.0525E-01
            -5.2572E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1630.34495759274        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  9.7566E-01  1.9037E-01  1.7233E+00  1.5878E+00  9.2630E-01  1.0572E+00  2.4972E+00  1.3368E+00  8.5801E-01  9.0904E-01
             1.0887E+00
 PARAMETER:  7.5356E-02 -1.5588E+00  6.4425E-01  5.6233E-01  2.3441E-02  1.5564E-01  1.0152E+00  3.9026E-01 -5.3134E-02  4.6308E-03
             1.8500E-01
 GRADIENT:  -7.3142E-01  3.0650E+00  1.9083E+00  1.6326E+01 -6.9416E+00 -2.4839E-01 -4.0657E-01  3.0014E-02 -8.0193E-01  8.9881E-01
            -3.0589E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1630.38820541707        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  9.7487E-01  1.3027E-01  1.7868E+00  1.6284E+00  9.2732E-01  1.0572E+00  2.9901E+00  1.3892E+00  8.4366E-01  9.1013E-01
             1.0894E+00
 PARAMETER:  7.4547E-02 -1.9382E+00  6.8044E-01  5.8763E-01  2.4542E-02  1.5560E-01  1.1953E+00  4.2875E-01 -7.0006E-02  5.8283E-03
             1.8559E-01
 GRADIENT:  -7.6139E-01  2.2749E+00  1.8662E+00  1.6194E+01 -6.9004E+00 -1.0928E-01 -2.6049E-01  1.9550E-01 -7.8071E-01  9.0064E-01
            -7.3670E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1630.48442294321        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1588
 NPARAMETR:  9.7418E-01  6.6953E-02  1.8609E+00  1.6692E+00  9.3094E-01  1.0571E+00  4.1781E+00  1.4478E+00  8.2869E-01  9.1042E-01
             1.0900E+00
 PARAMETER:  7.3838E-02 -2.6038E+00  7.2107E-01  6.1232E-01  2.8434E-02  1.5555E-01  1.5299E+00  4.7006E-01 -8.7914E-02  6.1548E-03
             1.8619E-01
 GRADIENT:  -5.1756E-01  1.1539E+00  1.3159E+00  1.1148E+01 -4.7447E+00  2.7578E-02 -1.0198E-01  2.7734E-01 -6.6337E-01  5.7335E-01
             1.1243E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1630.54968912702        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1763
 NPARAMETR:  9.7383E-01  3.0356E-02  1.9024E+00  1.6920E+00  9.3308E-01  1.0569E+00  6.3057E+00  1.4737E+00  8.2076E-01  9.1051E-01
             1.0902E+00
 PARAMETER:  7.3482E-02 -3.3948E+00  7.4314E-01  6.2590E-01  3.0741E-02  1.5535E-01  1.9415E+00  4.8780E-01 -9.7523E-02  6.2545E-03
             1.8639E-01
 GRADIENT:  -2.6529E-01  5.2716E-01  1.0130E+00  6.7573E+00 -2.6656E+00  3.9817E-02 -2.4178E-02 -1.7525E-03 -3.6997E-01  1.8075E-01
             4.5050E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1630.60696450538        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1940
 NPARAMETR:  9.7365E-01  1.0000E-02  1.9187E+00  1.7032E+00  9.3317E-01  1.0570E+00  1.1661E+01  1.4848E+00  8.1636E-01  9.1045E-01
             1.0902E+00
 PARAMETER:  7.3299E-02 -4.5555E+00  7.5163E-01  6.3253E-01  3.0829E-02  1.5548E-01  2.5562E+00  4.9531E-01 -1.0290E-01  6.1801E-03
             1.8640E-01
 GRADIENT:  -1.0992E-01  0.0000E+00  4.1869E-01  2.0815E+00 -9.1288E-01  1.6088E-01  1.3917E-02 -2.0798E-02 -2.6199E-01  3.4658E-02
             4.5644E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1630.73398207804        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2126             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7422E-01  1.0000E-02  1.9023E+00  1.6883E+00  9.3065E-01  1.0573E+00  1.1690E+01  1.4741E+00  8.1652E-01  9.0793E-01
             1.0896E+00
 PARAMETER:  7.3879E-02 -4.6001E+00  7.4308E-01  6.2374E-01  2.8130E-02  1.5576E-01  2.5588E+00  4.8802E-01 -1.0271E-01  3.4077E-03
             1.8583E-01
 GRADIENT:   3.6638E+02  0.0000E+00  8.1831E+00  1.0847E+03  7.7401E+00  6.7715E+01  7.8511E-01  1.2649E+00  1.6683E+01  3.9445E-01
             1.5159E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1630.73631857460        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2314             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7424E-01  1.0000E-02  1.8941E+00  1.6878E+00  9.2934E-01  1.0573E+00  1.1708E+01  1.4693E+00  8.1688E-01  9.0467E-01
             1.0894E+00
 PARAMETER:  7.3905E-02 -4.6001E+00  7.3874E-01  6.2343E-01  2.6720E-02  1.5573E-01  2.5603E+00  4.8478E-01 -1.0226E-01 -1.8870E-04
             1.8563E-01
 GRADIENT:   3.6649E+02  0.0000E+00  7.6823E+00  1.0838E+03  8.8954E+00  6.7678E+01  7.8986E-01  1.2241E+00  1.6700E+01  5.5862E-02
             1.3956E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1630.74171222412        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2501             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7423E-01  1.0000E-02  1.8908E+00  1.6876E+00  9.2750E-01  1.0573E+00  1.1720E+01  1.4654E+00  8.1699E-01  9.0604E-01
             1.0894E+00
 PARAMETER:  7.3887E-02 -4.6001E+00  7.3699E-01  6.2330E-01  2.4738E-02  1.5572E-01  2.5613E+00  4.8214E-01 -1.0212E-01  1.3258E-03
             1.8561E-01
 GRADIENT:   3.6647E+02  0.0000E+00  8.3649E+00  1.0836E+03  7.0965E+00  6.7674E+01  7.9131E-01  1.2417E+00  1.6706E+01  4.1490E-01
             1.4452E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1630.74331079124        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2689
 NPARAMETR:  9.7421E-01  1.0000E-02  1.8838E+00  1.6873E+00  9.2649E-01  1.0573E+00  1.1723E+01  1.4605E+00  8.1709E-01  9.0339E-01
             1.0893E+00
 PARAMETER:  7.3869E-02 -4.6001E+00  7.3328E-01  6.2312E-01  2.3646E-02  1.5569E-01  2.5616E+00  4.7875E-01 -1.0201E-01 -1.5993E-03
             1.8551E-01
 GRADIENT:   1.4557E+00  0.0000E+00 -2.7304E-01 -2.4435E+01  2.6244E+00  3.1445E-01  1.6222E-02  2.5823E-02 -1.9946E-02 -4.6523E-01
            -1.3065E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1630.74789099527        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2880             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7419E-01  1.0000E-02  1.8795E+00  1.6870E+00  9.2387E-01  1.0573E+00  1.1756E+01  1.4564E+00  8.1728E-01  9.0663E-01
             1.0895E+00
 PARAMETER:  7.3851E-02 -4.6001E+00  7.3099E-01  6.2295E-01  2.0819E-02  1.5570E-01  2.5643E+00  4.7596E-01 -1.0177E-01  1.9835E-03
             1.8570E-01
 GRADIENT:   3.6622E+02  0.0000E+00  8.8745E+00  1.0823E+03  5.3019E+00  6.7621E+01  7.9706E-01  1.3088E+00  1.6690E+01  8.1431E-01
             1.5722E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1630.74975041668        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3068
 NPARAMETR:  9.7418E-01  1.0000E-02  1.8756E+00  1.6868E+00  9.2312E-01  1.0573E+00  1.1756E+01  1.4531E+00  8.1734E-01  9.0547E-01
             1.0894E+00
 PARAMETER:  7.3844E-02 -4.6001E+00  7.2892E-01  6.2285E-01  2.0002E-02  1.5568E-01  2.5644E+00  4.7370E-01 -1.0170E-01  6.9500E-04
             1.8563E-01
 GRADIENT:   1.4724E+00  0.0000E+00  5.6351E-01 -2.4310E+01 -1.5356E-03  3.4215E-01  1.7425E-02  1.2875E-01  3.3271E-02  1.2645E-01
             3.1448E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1630.75194816021        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3260             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7417E-01  1.0000E-02  1.8690E+00  1.6865E+00  9.2218E-01  1.0572E+00  1.1763E+01  1.4480E+00  8.1744E-01  9.0326E-01
             1.0893E+00
 PARAMETER:  7.3832E-02 -4.6001E+00  7.2542E-01  6.2267E-01  1.8990E-02  1.5566E-01  2.5649E+00  4.7018E-01 -1.0157E-01 -1.7473E-03
             1.8555E-01
 GRADIENT:   3.6623E+02  0.0000E+00  8.3321E+00  1.0814E+03  6.6557E+00  6.7605E+01  7.9928E-01  1.1618E+00  1.6627E+01  4.2154E-01
             1.4382E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1630.75315094306        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3448
 NPARAMETR:  9.7416E-01  1.0000E-02  1.8662E+00  1.6864E+00  9.2137E-01  1.0572E+00  1.1771E+01  1.4455E+00  8.1752E-01  9.0326E-01
             1.0893E+00
 PARAMETER:  7.3824E-02 -4.6001E+00  7.2390E-01  6.2258E-01  1.8105E-02  1.5566E-01  2.5656E+00  4.6849E-01 -1.0148E-01 -1.7467E-03
             1.8555E-01
 GRADIENT:   1.4212E+00  0.0000E+00  1.9345E-01 -2.4349E+01  7.4958E-01  3.1809E-01  1.7213E-02  5.6032E-02  1.7332E-02 -8.4797E-02
            -3.2624E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1630.75763864315        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3641             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7412E-01  1.0000E-02  1.8436E+00  1.6853E+00  9.1666E-01  1.0572E+00  1.1803E+01  1.4291E+00  8.1797E-01  9.0403E-01
             1.0894E+00
 PARAMETER:  7.3775E-02 -4.6001E+00  7.1173E-01  6.2193E-01  1.2986E-02  1.5559E-01  2.5684E+00  4.5705E-01 -1.0093E-01 -8.9681E-04
             1.8563E-01
 GRADIENT:   3.6583E+02  0.0000E+00  7.7175E+00  1.0785E+03  6.5819E+00  6.7435E+01  8.0701E-01  1.3161E+00  1.6552E+01  8.5629E-01
             1.6239E+00

0ITERATION NO.:  104    OBJECTIVE VALUE:  -1630.75852213763        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     3780
 NPARAMETR:  9.7411E-01  1.0000E-02  1.8408E+00  1.6852E+00  9.1574E-01  1.0572E+00  1.1844E+01  1.4314E+00  8.1805E-01  9.0250E-01
             1.0899E+00
 PARAMETER:  7.3774E-02 -4.6001E+00  7.1250E-01  6.2193E-01  1.2451E-02  1.5560E-01  2.5686E+00  4.5457E-01 -1.0090E-01 -3.5873E-03
             1.8543E-01
 GRADIENT:   1.2695E-03  0.0000E+00  2.9415E-01  7.4690E-02  2.1584E-01 -5.3950E-03 -4.2071E-04 -1.0482E-01 -8.0173E-03 -3.8428E-02
            -1.0090E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3780
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7514E-06 -2.1482E-04 -2.7425E-02 -5.6515E-03 -3.8494E-02
 SE:             2.9808E-02  1.7565E-03  1.7952E-02  2.9268E-02  1.9849E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9993E-01  9.0267E-01  1.2660E-01  8.4689E-01  5.2462E-02

 ETASHRINKSD(%)  1.4036E-01  9.4115E+01  3.9857E+01  1.9471E+00  3.3502E+01
 ETASHRINKVR(%)  2.8053E-01  9.9654E+01  6.3828E+01  3.8563E+00  5.5780E+01
 EBVSHRINKSD(%)  4.5303E-01  9.4234E+01  4.3032E+01  2.2995E+00  3.1148E+01
 EBVSHRINKVR(%)  9.0400E-01  9.9668E+01  6.7547E+01  4.5461E+00  5.2595E+01
 RELATIVEINF(%)  9.4836E+01  6.9334E-03  6.7453E+00  2.4228E+00  5.5604E+00
 EPSSHRINKSD(%)  4.3969E+01
 EPSSHRINKVR(%)  6.8605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1630.7585221376296     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -895.60769557389142     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1630.759       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.00E-02  1.85E+00  1.69E+00  9.16E-01  1.06E+00  1.18E+01  1.43E+00  8.18E-01  9.02E-01  1.09E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.41E+00  0.00E+00  5.50E+01
 
 TH 4
+       -7.70E+00  0.00E+00 -1.19E+01  5.59E+02
 
 TH 5
+       -1.64E+00  0.00E+00 -1.67E+02 -5.57E+01  7.93E+02
 
 TH 6
+        7.91E-01  0.00E+00  1.30E-01 -2.37E+00 -2.00E+00  1.74E+02
 
 TH 7
+        9.54E-04  0.00E+00  2.01E-03 -3.26E-03 -2.21E-03 -3.22E-04  1.46E-03
 
 TH 8
+        1.96E-01  0.00E+00 -1.45E+01 -2.44E+00 -5.68E+00 -5.86E-02  1.97E-03  1.84E+01
 
 TH 9
+        2.50E+00  0.00E+00  5.14E+00 -7.09E-01 -1.32E+00 -4.30E-01  3.78E-02  9.28E-01  2.73E+02
 
 TH10
+       -3.20E-01  0.00E+00  7.00E-01 -6.76E-01 -8.80E+01  7.43E-01  6.85E-03  1.52E+01  1.54E+00  6.89E+01
 
 TH11
+       -8.23E+00  0.00E+00 -7.19E+00 -6.69E+00 -9.55E+00  2.43E+00  1.71E-03  9.95E+00  5.66E+00  1.34E+01  1.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       61.906
Stop Time:
Wed Sep 29 11:07:45 CDT 2021
