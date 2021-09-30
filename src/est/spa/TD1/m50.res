Wed Sep 29 18:13:52 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1689.67915742881        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2349E+02 -1.0712E+01 -2.3406E+01  8.0728E+00  4.6391E-01  3.7702E+01  9.8741E+00  1.8242E+01  1.1047E+01  2.3526E+01
             7.2089E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1693.74735348736        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9235E-01  1.0366E+00  1.2004E+00  1.0210E+00  1.0778E+00  1.0254E+00  9.3769E-01  8.5970E-01  9.5970E-01  8.6147E-01
             9.8689E-01
 PARAMETER:  9.2316E-02  1.3599E-01  2.8263E-01  1.2075E-01  1.7490E-01  1.2505E-01  3.5664E-02 -5.1174E-02  5.8865E-02 -4.9111E-02
             8.6808E-02
 GRADIENT:  -7.7386E+00  1.8013E+01  2.7627E+01 -1.7386E+01 -1.5806E+01  5.8778E+00  2.1015E+00  1.8689E+00 -1.0940E+01 -2.0305E+01
            -1.5343E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.82931577573        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.9414E-01  9.0922E-01  1.2709E+00  1.1165E+00  1.1001E+00  9.7104E-01  6.6838E-01  3.8626E-01  9.7284E-01  1.1315E+00
             1.0451E+00
 PARAMETER:  9.4118E-02  4.8361E-03  3.3971E-01  2.1018E-01  1.9541E-01  7.0617E-02 -3.0290E-01 -8.5126E-01  7.2464E-02  2.2353E-01
             1.4407E-01
 GRADIENT:  -3.1141E+00  1.7061E+01 -4.8765E+00  1.9811E+01  6.1690E+00 -1.4940E+01 -1.5455E+00  4.6877E-01 -1.1503E+01  7.1092E+00
             1.6209E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1699.94483346744        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.9443E-01  7.1698E-01  1.1567E+00  1.2308E+00  9.4890E-01  1.0103E+00  9.1349E-01  2.8016E-01  9.0562E-01  9.8792E-01
             9.9850E-01
 PARAMETER:  9.4410E-02 -2.3271E-01  2.4560E-01  3.0765E-01  4.7552E-02  1.1027E-01  9.5211E-03 -1.1724E+00  8.6602E-04  8.7843E-02
             9.8494E-02
 GRADIENT:   6.8891E-01  1.3276E+01  6.9004E+00  1.7121E+01 -1.3825E+01  1.3151E+00 -2.1204E-01  1.2783E-01 -5.2963E-01 -8.3902E-01
            -1.4207E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.74412436321        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.9088E-01  4.9002E-01  1.2423E+00  1.3718E+00  9.1588E-01  9.9704E-01  9.0009E-01  1.2107E-01  8.3603E-01  1.0376E+00
             9.9881E-01
 PARAMETER:  9.0834E-02 -6.1332E-01  3.1694E-01  4.1614E-01  1.2131E-02  9.7035E-02 -5.2555E-03 -2.0114E+00 -7.9092E-02  1.3692E-01
             9.8806E-02
 GRADIENT:  -1.7596E-01  8.6635E+00  6.5558E+00  1.8084E+01 -1.1613E+01 -2.4336E+00  2.1596E-02 -5.4538E-03 -1.9330E+00  1.0607E+00
            -1.6936E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.21331932702        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.8784E-01  2.8738E-01  1.3003E+00  1.4904E+00  8.8892E-01  9.9767E-01  6.9858E-01  2.7287E-02  7.8642E-01  1.0517E+00
             1.0007E+00
 PARAMETER:  8.7766E-02 -1.1470E+00  3.6261E-01  4.9904E-01 -1.7744E-02  9.7669E-02 -2.5871E-01 -3.5013E+00 -1.4026E-01  1.5038E-01
             1.0071E-01
 GRADIENT:   4.5118E-01  1.9721E+00  1.0514E+00  3.8340E+00  1.0466E+00 -7.7508E-01  2.3993E-02 -9.3407E-04  3.9296E-01 -4.2259E-01
            -1.0239E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.22828261135        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.8684E-01  2.3320E-01  1.2944E+00  1.5256E+00  8.6981E-01  9.9750E-01  6.3214E-01  1.4314E-02  7.6789E-01  1.0457E+00
             1.0006E+00
 PARAMETER:  8.6751E-02 -1.3558E+00  3.5806E-01  5.2237E-01 -3.9481E-02  9.7502E-02 -3.5864E-01 -4.1465E+00 -1.6411E-01  1.4470E-01
             1.0064E-01
 GRADIENT:   1.1643E-01  2.6125E+00  1.2738E+00  1.0984E+01 -1.1819E+00 -4.8263E-01 -4.0895E-03 -2.3344E-04 -1.0982E+00 -3.4534E-01
            -9.7678E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.32005204466        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1248
 NPARAMETR:  9.8721E-01  2.2383E-01  1.2934E+00  1.5211E+00  8.6706E-01  9.9872E-01  4.2964E-01  2.1935E-01  7.6791E-01  1.0474E+00
             1.0016E+00
 PARAMETER:  8.7125E-02 -1.3969E+00  3.5728E-01  5.1946E-01 -4.2644E-02  9.8723E-02 -7.4482E-01 -1.4171E+00 -1.6409E-01  1.4634E-01
             1.0165E-01
 GRADIENT:   1.4838E+00  2.8380E-01 -5.7791E-02 -1.2256E+01 -1.6066E-02  1.1333E-01  3.8908E-03  2.2191E-02 -2.8618E-01  2.0331E+00
             1.1736E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.33042531044        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.8464E-01  2.1542E-01  1.2915E+00  1.5259E+00  8.6235E-01  9.9738E-01  1.8009E-01  2.5063E-01  7.6662E-01  1.0332E+00
             9.9834E-01
 PARAMETER:  8.4524E-02 -1.4352E+00  3.5580E-01  5.2256E-01 -4.8096E-02  9.7375E-02 -1.6143E+00 -1.2838E+00 -1.6576E-01  1.3266E-01
             9.8341E-02
 GRADIENT:  -3.9116E+00  4.4775E-01  1.5040E+00 -1.2147E+01 -1.2676E+00 -3.9662E-01  6.3050E-04 -1.2772E-02 -1.9091E-01 -1.3041E-02
            -7.3824E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.33049452061        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1606
 NPARAMETR:  9.8430E-01  2.1264E-01  1.2901E+00  1.5270E+00  8.6095E-01  9.9718E-01  1.2938E-01  2.5264E-01  7.6590E-01  1.0317E+00
             9.9826E-01
 PARAMETER:  8.4174E-02 -1.4481E+00  3.5473E-01  5.2334E-01 -4.9724E-02  9.7175E-02 -1.9450E+00 -1.2758E+00 -1.6670E-01  1.3119E-01
             9.8257E-02
 GRADIENT:  -4.5785E+00  3.6987E-01  1.4721E+00 -1.3046E+01 -1.1970E+00 -4.5900E-01  6.0961E-04 -1.2817E-02 -1.6807E-01 -1.2384E-01
            -7.5278E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1701.33053529008        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1786
 NPARAMETR:  9.8402E-01  2.0999E-01  1.2883E+00  1.5282E+00  8.5948E-01  9.9701E-01  8.8076E-02  2.4938E-01  7.6510E-01  1.0309E+00
             9.9824E-01
 PARAMETER:  8.3890E-02 -1.4607E+00  3.5335E-01  5.2408E-01 -5.1433E-02  9.7002E-02 -2.3296E+00 -1.2888E+00 -1.6775E-01  1.3045E-01
             9.8236E-02
 GRADIENT:  -5.1074E+00  3.0654E-01  1.4188E+00 -1.3813E+01 -1.1374E+00 -5.0959E-01  4.6631E-04 -1.2672E-02 -1.8674E-01 -1.6267E-01
            -7.4322E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1701.33054878109        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  9.8385E-01  2.0822E-01  1.2870E+00  1.5289E+00  8.5842E-01  9.9690E-01  6.5885E-02  2.4584E-01  7.6456E-01  1.0305E+00
             9.9824E-01
 PARAMETER:  8.3720E-02 -1.4692E+00  3.5230E-01  5.2458E-01 -5.2656E-02  9.6896E-02 -2.6198E+00 -1.3031E+00 -1.6846E-01  1.3008E-01
             9.8235E-02
 GRADIENT:  -5.4178E+00  2.6854E-01  1.3835E+00 -1.4281E+01 -1.1131E+00 -5.3971E-01  3.5017E-04 -1.2513E-02 -1.9394E-01 -1.7278E-01
            -7.3137E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1701.33055525281        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2151
 NPARAMETR:  9.8373E-01  2.0686E-01  1.2859E+00  1.5295E+00  8.5759E-01  9.9683E-01  5.1553E-02  2.4255E-01  7.6414E-01  1.0303E+00
             9.9824E-01
 PARAMETER:  8.3600E-02 -1.4757E+00  3.5142E-01  5.2497E-01 -5.3634E-02  9.6821E-02 -2.8651E+00 -1.3165E+00 -1.6901E-01  1.2983E-01
             9.8237E-02
 GRADIENT:  -5.6333E+00  2.4198E-01  1.3540E+00 -1.4613E+01 -1.0967E+00 -5.6027E-01  2.6463E-04 -1.2398E-02 -2.0100E-01 -1.7703E-01
            -7.2177E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1701.34507393930        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  9.8675E-01  2.0573E-01  1.2847E+00  1.5314E+00  8.5773E-01  9.9859E-01  3.6391E-02  2.4488E-01  7.6448E-01  1.0299E+00
             9.9987E-01
 PARAMETER:  8.6663E-02 -1.4812E+00  3.5054E-01  5.2618E-01 -5.3471E-02  9.8591E-02 -3.2134E+00 -1.3070E+00 -1.6856E-01  1.2942E-01
             9.9872E-02
 GRADIENT:   1.1535E+00  2.8261E-01 -9.0759E-02 -1.1906E+01  6.4274E-01  1.6048E-01  1.7428E-04 -2.8031E-03  1.5257E-01 -1.3812E-01
             6.8882E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1701.34588193010        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2501
 NPARAMETR:  9.8680E-01  2.0384E-01  1.2851E+00  1.5324E+00  8.5729E-01  9.9841E-01  3.3089E-02  2.5291E-01  7.6394E-01  1.0305E+00
             9.9972E-01
 PARAMETER:  8.6716E-02 -1.4904E+00  3.5087E-01  5.2681E-01 -5.3974E-02  9.8409E-02 -3.3086E+00 -1.2747E+00 -1.6927E-01  1.3000E-01
             9.9718E-02
 GRADIENT:   1.3577E+00  2.5007E-01 -9.4912E-02 -1.2418E+01  4.2281E-01  1.0401E-01  1.5502E-04  7.1324E-03  1.8542E-01  1.4140E-01
             1.2776E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1701.34699243146        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2685
 NPARAMETR:  9.8675E-01  2.0075E-01  1.2860E+00  1.5342E+00  8.5672E-01  9.9836E-01  2.8687E-02  2.5037E-01  7.6299E-01  1.0303E+00
             9.9959E-01
 PARAMETER:  8.6660E-02 -1.5057E+00  3.5151E-01  5.2800E-01 -5.4645E-02  9.8360E-02 -3.4513E+00 -1.2848E+00 -1.7051E-01  1.2982E-01
             9.9586E-02
 GRADIENT:   1.3667E+00  2.3990E-01  1.2196E-01 -1.2683E+01  2.8424E-01  1.0708E-01  1.2612E-04  1.1161E-03  2.0576E-01  3.6571E-02
             5.7956E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1701.34716945455        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2864
 NPARAMETR:  9.8650E-01  1.9683E-01  1.2869E+00  1.5377E+00  8.5626E-01  9.9823E-01  2.2580E-02  2.5127E-01  7.6204E-01  1.0304E+00
             9.9961E-01
 PARAMETER:  8.6405E-02 -1.5254E+00  3.5223E-01  5.3026E-01 -5.5183E-02  9.8232E-02 -3.6907E+00 -1.2812E+00 -1.7176E-01  1.2990E-01
             9.9607E-02
 GRADIENT:   9.4622E-01  3.9400E-01 -2.0147E-01 -1.0332E+01  5.9072E-01  7.8828E-02  8.7835E-05  1.5501E-04  3.4435E-01  1.7570E-02
             7.6148E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1701.34852376524        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3049
 NPARAMETR:  9.8656E-01  1.9470E-01  1.2873E+00  1.5383E+00  8.5582E-01  9.9823E-01  2.0282E-02  2.5279E-01  7.6123E-01  1.0302E+00
             9.9959E-01
 PARAMETER:  8.6469E-02 -1.5363E+00  3.5258E-01  5.3068E-01 -5.5696E-02  9.8228E-02 -3.7980E+00 -1.2752E+00 -1.7282E-01  1.2979E-01
             9.9591E-02
 GRADIENT:   1.1891E+00  2.7595E-01 -7.0830E-02 -1.1959E+01  5.5332E-01  9.5925E-02  7.4618E-05  5.8431E-04  2.8511E-01  2.9573E-02
             1.2852E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1701.34928070910        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3235
 NPARAMETR:  9.8659E-01  1.9246E-01  1.2879E+00  1.5394E+00  8.5529E-01  9.9820E-01  1.7461E-02  2.5467E-01  7.6027E-01  1.0301E+00
             9.9957E-01
 PARAMETER:  8.6497E-02 -1.5479E+00  3.5301E-01  5.3141E-01 -5.6319E-02  9.8201E-02 -3.9478E+00 -1.2678E+00 -1.7408E-01  1.2967E-01
             9.9568E-02
 GRADIENT:   1.3478E+00  2.5237E-01  1.3344E-01 -1.2605E+01  2.9898E-01  1.0311E-01  5.8695E-05  6.0988E-04  1.6857E-01  5.2733E-02
            -6.7670E-04

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1701.34976864979        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3426
 NPARAMETR:  9.8662E-01  1.9142E-01  1.2874E+00  1.5397E+00  8.5468E-01  9.9820E-01  2.1269E-02  2.5632E-01  7.5982E-01  1.0294E+00
             9.9953E-01
 PARAMETER:  8.6534E-02 -1.5533E+00  3.5260E-01  5.3160E-01 -5.7028E-02  9.8195E-02 -3.7505E+00 -1.2613E+00 -1.7467E-01  1.2900E-01
             9.9534E-02
 GRADIENT:   1.4753E+00  2.1246E-01  2.1367E-01 -1.3283E+01  2.0534E-01  1.0918E-01  7.7525E-05  1.1299E-03  1.0899E-01  2.2881E-02
            -8.4904E-03

0ITERATION NO.:   98    OBJECTIVE VALUE:  -1701.34981144284        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     3522
 NPARAMETR:  9.8654E-01  1.9100E-01  1.2874E+00  1.5405E+00  8.5465E-01  9.9816E-01  1.9725E-02  2.5623E-01  7.5978E-01  1.0294E+00
             9.9954E-01
 PARAMETER:  8.6451E-02 -1.5555E+00  3.5260E-01  5.3208E-01 -5.7061E-02  9.8162E-02 -3.8259E+00 -1.2617E+00 -1.7473E-01  1.2900E-01
             9.9542E-02
 GRADIENT:  -1.5305E-01  7.7552E-02  7.8125E-02  6.6568E-01  3.3711E-01 -8.7032E-03 -1.0766E-06 -8.9071E-04  9.3969E-02  1.0215E-02
            -7.3969E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3522
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1036E-05 -1.4079E-04 -5.4491E-03 -5.9170E-03 -2.3923E-02
 SE:             2.9830E-02  7.1578E-05  4.4368E-03  2.9200E-02  2.4962E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9863E-01  4.9185E-02  2.1939E-01  8.3942E-01  3.3786E-01

 ETASHRINKSD(%)  6.5558E-02  9.9760E+01  8.5136E+01  2.1761E+00  1.6375E+01
 ETASHRINKVR(%)  1.3107E-01  9.9999E+01  9.7791E+01  4.3049E+00  3.0069E+01
 EBVSHRINKSD(%)  4.0926E-01  9.9778E+01  8.5309E+01  2.2925E+00  1.3464E+01
 EBVSHRINKVR(%)  8.1684E-01  1.0000E+02  9.7842E+01  4.5324E+00  2.5116E+01
 RELATIVEINF(%)  9.7215E+01  2.4304E-05  2.6779E-01  6.1407E+00  6.0707E+00
 EPSSHRINKSD(%)  4.1750E+01
 EPSSHRINKVR(%)  6.6069E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.3498114428403     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.19898487910211     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.30
 Elapsed covariance  time in seconds:     5.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.350       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  1.91E-01  1.29E+00  1.54E+00  8.55E-01  9.98E-01  1.97E-02  2.56E-01  7.60E-01  1.03E+00  1.00E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.03E-02  2.55E-01  2.24E-01  1.56E-01  1.14E-01  7.27E-02  1.50E-01  8.21E-01  1.01E-01  1.37E-01  6.78E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.18E-04
 
 TH 2
+        1.61E-03  6.52E-02
 
 TH 3
+       -3.40E-04 -8.59E-03  5.03E-02
 
 TH 4
+       -9.34E-04 -3.84E-02  1.05E-02  2.44E-02
 
 TH 5
+        3.62E-04  1.28E-02  1.97E-02 -5.07E-03  1.30E-02
 
 TH 6
+       -2.43E-04  1.47E-03  2.96E-05 -7.03E-04  4.51E-04  5.28E-03
 
 TH 7
+       -2.30E-04 -1.37E-02  1.03E-02  8.60E-03 -1.35E-03 -1.13E-03  2.24E-02
 
 TH 8
+       -2.00E-03 -5.72E-02  5.94E-02  3.71E-02  7.72E-04 -1.05E-02  1.20E-01  6.73E-01
 
 TH 9
+        4.26E-04  1.93E-02 -2.74E-05 -1.11E-02  4.71E-03 -4.56E-04 -4.09E-03 -1.54E-02  1.02E-02
 
 TH10
+       -2.69E-04  6.53E-03  1.23E-02 -2.34E-03  9.24E-03  1.98E-03 -1.08E-02 -4.90E-02  2.04E-03  1.87E-02
 
 TH11
+        1.69E-05  1.42E-03  2.13E-03 -2.74E-04  1.56E-03  3.82E-04 -2.69E-03 -1.91E-02  4.07E-04  9.75E-04  4.59E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.03E-02
 
 TH 2
+        2.08E-01  2.55E-01
 
 TH 3
+       -5.00E-02 -1.50E-01  2.24E-01
 
 TH 4
+       -1.97E-01 -9.61E-01  3.00E-01  1.56E-01
 
 TH 5
+        1.05E-01  4.42E-01  7.71E-01 -2.85E-01  1.14E-01
 
 TH 6
+       -1.11E-01  7.94E-02  1.81E-03 -6.19E-02  5.45E-02  7.27E-02
 
 TH 7
+       -5.08E-02 -3.58E-01  3.07E-01  3.68E-01 -7.93E-02 -1.04E-01  1.50E-01
 
 TH 8
+       -8.04E-02 -2.73E-01  3.23E-01  2.90E-01  8.26E-03 -1.75E-01  9.76E-01  8.21E-01
 
 TH 9
+        1.40E-01  7.51E-01 -1.21E-03 -7.06E-01  4.10E-01 -6.23E-02 -2.72E-01 -1.86E-01  1.01E-01
 
 TH10
+       -6.49E-02  1.87E-01  4.00E-01 -1.09E-01  5.93E-01  2.00E-01 -5.26E-01 -4.36E-01  1.48E-01  1.37E-01
 
 TH11
+        8.24E-03  8.22E-02  1.40E-01 -2.59E-02  2.03E-01  7.76E-02 -2.65E-01 -3.44E-01  5.96E-02  1.05E-01  6.78E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.93E+03
 
 TH 2
+        8.50E+01  4.30E+02
 
 TH 3
+        5.99E+02  1.68E+02  6.82E+02
 
 TH 4
+       -1.18E+02  4.88E+02 -1.43E+02  8.47E+02
 
 TH 5
+       -1.25E+03 -4.97E+02 -1.33E+03  6.96E+01  2.89E+03
 
 TH 6
+        4.73E+02  6.20E+00  3.03E+02 -8.73E+01 -5.60E+02  4.68E+02
 
 TH 7
+       -2.49E+03 -2.99E+02 -1.94E+03  3.81E+02  3.70E+03 -1.51E+03  9.00E+03
 
 TH 8
+        4.12E+02  4.77E+01  3.02E+02 -5.95E+01 -5.95E+02  2.46E+02 -1.45E+03  2.39E+02
 
 TH 9
+       -1.22E+02 -1.08E+02 -1.59E+02  1.99E+01  2.69E+02 -4.75E+01  5.53E+02 -8.50E+01  2.82E+02
 
 TH10
+       -2.01E+02  1.10E+01 -2.10E+02  6.83E+01  2.35E+02 -2.08E+02  1.09E+03 -1.61E+02  8.19E+01  3.00E+02
 
 TH11
+        3.80E+02  1.76E+01  2.31E+02 -9.12E+01 -5.45E+02  1.91E+02 -1.18E+03  2.07E+02 -5.18E+01 -6.65E+01  4.61E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       49.889
Stop Time:
Wed Sep 29 18:14:44 CDT 2021
