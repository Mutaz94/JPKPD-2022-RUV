Sat Oct 23 15:46:49 CDT 2021
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
$DATA ../../../../data/SD1/TD1/dat38.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3389.33139410697        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6134E+02  6.7161E+01  1.5022E+02  6.7602E+01  6.8890E+01  1.5620E+01 -6.0216E+01 -4.6640E+02 -1.4095E+02 -1.0851E+01
            -2.0024E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3590.00321702476        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      151
 NPARAMETR:  9.9189E-01  1.0349E+00  9.8023E-01  9.5445E-01  1.0429E+00  8.3730E-01  1.1639E+00  2.0264E+00  1.0551E+00  9.7446E-01
             1.2517E+00
 PARAMETER:  9.1859E-02  1.3434E-01  8.0031E-02  5.3385E-02  1.4197E-01 -7.7575E-02  2.5180E-01  8.0627E-01  1.5360E-01  7.4133E-02
             3.2452E-01
 GRADIENT:  -3.6995E+00 -5.3426E+01 -3.0238E+00 -4.1225E+01 -4.5057E+01 -1.2784E+02 -2.1684E+01 -1.0510E+02 -1.6719E+01  1.1246E+01
             2.9672E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3608.67047711478        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  9.6869E-01  1.0958E+00  1.3170E+00  8.8773E-01  1.3674E+00  8.3852E-01  1.3913E+00  2.6967E+00  7.8676E-01  1.4472E+00
             1.2103E+00
 PARAMETER:  6.8189E-02  1.9151E-01  3.7532E-01 -1.9086E-02  4.1290E-01 -7.6121E-02  4.3023E-01  1.0920E+00 -1.3984E-01  4.6967E-01
             2.9087E-01
 GRADIENT:  -7.7379E+01 -1.2054E+02  4.4972E+00 -3.0742E+01  4.6836E+01 -1.2953E+02 -1.0718E+01 -7.3793E+01 -1.4587E+01  4.4484E+01
             2.3349E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3632.64174608864        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.5382E-01  1.2214E+00  1.1837E+00  8.4372E-01  1.3789E+00  8.6620E-01  1.2826E+00  2.9168E+00  7.4953E-01  1.3885E+00
             1.1654E+00
 PARAMETER:  5.2717E-02  3.0000E-01  2.6862E-01 -6.9937E-02  4.2132E-01 -4.3644E-02  3.4889E-01  1.1705E+00 -1.8831E-01  4.2824E-01
             2.5310E-01
 GRADIENT:  -1.1723E+02 -7.5339E+01 -7.5723E-01 -1.6772E+01  3.2046E+01 -1.1305E+02 -1.6125E+01 -4.6417E+01 -2.2234E+01  2.2656E+01
             1.6777E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3645.04778102721        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      652
 NPARAMETR:  9.5441E-01  1.2227E+00  1.1999E+00  8.4343E-01  1.3810E+00  1.0849E+00  1.2843E+00  2.9306E+00  7.4901E-01  1.3906E+00
             1.1643E+00
 PARAMETER:  5.3341E-02  3.0110E-01  2.8226E-01 -7.0281E-02  4.2284E-01  1.8144E-01  3.5018E-01  1.1752E+00 -1.8900E-01  4.2977E-01
             2.5215E-01
 GRADIENT:  -7.2609E+01 -7.4874E+01  1.4821E+00 -1.8096E+01  3.1034E+01 -7.5110E-01 -1.6022E+01 -4.5972E+01 -2.2211E+01  2.3018E+01
             1.6746E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3649.38104052127        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0239E+00  1.2271E+00  1.2486E+00  8.5065E-01  1.3758E+00  1.0978E+00  1.2866E+00  3.1428E+00  7.4954E-01  1.3865E+00
             1.1575E+00
 PARAMETER:  1.2359E-01  3.0462E-01  3.2199E-01 -6.1754E-02  4.1904E-01  1.9330E-01  3.5201E-01  1.2451E+00 -1.8829E-01  4.2680E-01
             2.4626E-01
 GRADIENT:   5.8780E+01 -6.4242E+01 -1.4954E+00 -7.0692E+00  1.9617E+01  4.9429E+00 -1.6797E+01 -3.0820E+01 -2.1119E+01  2.8427E+01
             1.5560E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3655.43613560985        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1011             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0117E+00  1.2419E+00  1.2001E+00  8.5057E-01  1.3765E+00  1.0791E+00  1.3015E+00  3.1818E+00  8.0435E-01  1.3578E+00
             1.1373E+00
 PARAMETER:  1.1166E-01  3.1666E-01  2.8239E-01 -6.1843E-02  4.1952E-01  1.7617E-01  3.6351E-01  1.2574E+00 -1.1773E-01  4.0587E-01
             2.2865E-01
 GRADIENT:   4.3863E+02  2.3418E+02  7.3901E+00  5.8380E+01  1.6226E+02  7.7357E+01  4.4811E+01  2.2329E+01 -1.2331E+01  4.0000E+01
             1.3358E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3665.59356288464        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1104
 NPARAMETR:  9.9543E-01  1.2443E+00  1.2284E+00  8.5579E-01  1.3785E+00  1.0683E+00  1.2920E+00  3.3650E+00  9.9010E-01  1.2261E+00
             1.0618E+00
 PARAMETER:  9.5420E-02  3.1860E-01  3.0569E-01 -5.5725E-02  4.2096E-01  1.6607E-01  3.5619E-01  1.3134E+00  9.0047E-02  3.0387E-01
             1.5993E-01
 GRADIENT:   8.4197E+00 -3.8424E+01 -8.7158E+00 -5.7378E+00  4.1599E+01 -4.9012E+00  4.9120E+00 -8.7965E+00 -8.6699E-01  1.2944E+01
             5.6410E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3666.43162365774        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1269
 NPARAMETR:  9.9524E-01  1.2446E+00  1.2877E+00  8.5577E-01  1.3775E+00  1.0621E+00  1.2466E+00  3.3589E+00  1.0052E+00  1.1506E+00
             1.0619E+00
 PARAMETER:  9.5231E-02  3.1878E-01  3.5284E-01 -5.5755E-02  4.2029E-01  1.6028E-01  3.2039E-01  1.3116E+00  1.0518E-01  2.4025E-01
             1.6007E-01
 GRADIENT:   8.0406E+00 -3.7303E+01 -2.8183E+00 -1.2206E+01  4.7902E+01 -7.2606E+00 -1.5812E+00 -9.9263E+00 -2.9998E+00  2.6223E+00
            -2.1624E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3667.29608210001        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1458             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9327E-01  1.2436E+00  1.2880E+00  8.5557E-01  1.3484E+00  1.0785E+00  1.2450E+00  3.3485E+00  1.0200E+00  1.1495E+00
             1.0623E+00
 PARAMETER:  9.3243E-02  3.1803E-01  3.5311E-01 -5.5983E-02  3.9892E-01  1.7558E-01  3.1915E-01  1.3085E+00  1.1982E-01  2.3930E-01
             1.6045E-01
 GRADIENT:   4.1512E+02  2.9486E+02  1.2289E+01  4.7413E+01  1.7137E+02  9.2987E+01  4.4037E+01  3.9125E+01  4.8993E+00  1.6938E+01
             3.6767E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3667.44639455352        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1625
 NPARAMETR:  9.9328E-01  1.2436E+00  1.2881E+00  8.5728E-01  1.3411E+00  1.0718E+00  1.2446E+00  3.3477E+00  1.0217E+00  1.1495E+00
             1.0623E+00
 PARAMETER:  9.3260E-02  3.1801E-01  3.5314E-01 -5.3988E-02  3.9351E-01  1.6939E-01  3.1878E-01  1.3083E+00  1.2147E-01  2.3934E-01
             1.6046E-01
 GRADIENT:   4.2342E+00 -2.8434E+01 -4.6702E-01 -1.9745E+01  2.0254E+01 -3.4838E+00 -5.6756E-01 -8.2521E+00 -1.9274E+00  6.9261E+00
             1.4233E+00

0ITERATION NO.:   51    OBJECTIVE VALUE:  -3667.44639455352        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1654
 NPARAMETR:  9.9333E-01  1.2434E+00  1.2879E+00  8.5725E-01  1.3411E+00  1.0737E+00  1.2452E+00  3.3472E+00  1.0217E+00  1.1495E+00
             1.0624E+00
 PARAMETER:  9.3260E-02  3.1801E-01  3.5314E-01 -5.3988E-02  3.9351E-01  1.6939E-01  3.1878E-01  1.3083E+00  1.2147E-01  2.3934E-01
             1.6046E-01
 GRADIENT:  -1.7545E+04  5.4852E+03  1.9936E+02  1.6827E+04  2.5599E+01 -3.8120E+00 -5.3925E-01  1.1932E+03 -1.8038E+00 -3.5052E-01
            -4.3095E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         3.4         1.0         1.8         2.9         5.1         5.5
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1654
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.4894E-04 -2.0435E-02 -2.6149E-02  4.1818E-02 -6.6985E-02
 SE:             3.0207E-02  2.4880E-02  2.5271E-02  2.4576E-02  2.2086E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7494E-01  4.1145E-01  3.0079E-01  8.8841E-02  2.4219E-03

 ETASHRINKSD(%)  1.0000E-10  1.6650E+01  1.5339E+01  1.7666E+01  2.6010E+01
 ETASHRINKVR(%)  1.0000E-10  3.0527E+01  2.8325E+01  3.2212E+01  4.5255E+01
 EBVSHRINKSD(%)  2.6356E-01  1.7836E+01  1.7311E+01  2.1589E+01  2.1307E+01
 EBVSHRINKVR(%)  5.2643E-01  3.2491E+01  3.1626E+01  3.8517E+01  3.8074E+01
 RELATIVEINF(%)  9.9472E+01  3.5461E+01  6.1489E+01  3.2963E+01  3.9116E+01
 EPSSHRINKSD(%)  2.2979E+01
 EPSSHRINKVR(%)  4.0678E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3667.4463945535167     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2013.3570347851060     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3667.446       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.24E+00  1.29E+00  8.57E-01  1.34E+00  1.07E+00  1.24E+00  3.35E+00  1.02E+00  1.15E+00  1.06E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,      141.814
Stop Time:
Sat Oct 23 15:47:11 CDT 2021
