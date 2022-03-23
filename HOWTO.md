# Using emDNA
The emDNA software package is developed as a tool to conduct elastic energy minimization calculations on a collection of DNA base pairs at the base-pair step level. 
Note: Arguments, flags and options are passed as with any standard UNIX programs, that is:
```
$ emDNA_tool --argument=value --flag --optional-parameter=some_value --optional-flag
```

## Data formats
emDNA can read and write data base-pair collection data in three formats, each with their own benefits:
- x3DNA base-pair step parameters format ([x3DNAparams](#x3DNAparams))<br>
- x3DNA base-pair reference frames format ([x3DNAbp](#x3DNAbp))<br>
- base pair list format ([bplist](#bplist))<br>

### x3DNAparams
The x3DNAparams format corresponds to a list of all the step parameters. The file usually contains a header indicating the type of the parameters. emDNA only uses the rigid-body base-pair parameters and hence the base parameters are not relevant; therefore, it is important that the value in the second line of the header is set to 0. Also, the header is not required for parsing.

An example of the parameter format 
```
#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist
T-A     -0.083    -0.197    -0.212     1.146   -12.278    -5.051    -0.264    -0.347     3.194     2.012    -1.575    35.719
```
The parsing of such file in emDNA expects that every line which is not describing a set of step parameters is either empty or starting with a #.
Details surrounding the parameters file can be found at [w3DNA 2.0](http://web.x3dna.org/index.php/rebuild)

### x3DNAbp
The x3DNAbp format corresponds to a description of every base pair in the collection as a series of origin values and reference frames for each base pair. 

An example of a base pair format:
```
...     5 A-T
   -5.3173    18.3644   -13.2982  # origin
   -0.0353    -0.1103     0.9933  # d1-axis
   -0.0890    -0.9896    -0.1131  # d2-axis
    0.9954    -0.0923     0.0251  # d3-axis
```
The first line indicates the index of the base pair (numbering starts at 1). The first line always starts with the ... characters and always contains the base-pair sequence after the index (that is, the two nucleotides forming the base pair). 
The second line describe the coordinates of the base-pair origin with respect to a global Cartesian coordinate system.
The remaining three lines describe the coordinate frame of the base pair. While the data after the # is not relevant and can be omitted, the third line describes the d1 or "short" axis that points into the major groove; the fourth describes the d2 ("long" axis) that points toward the backbone chain of the coding strand, and the fifth describe the d3 axis, which is cross product of d1 and d2 and points in the 5\`-3\` direction.

The parsing of such file in emDNA expects that every line which is not describing a base pair is either empty or starting with a #.

### bplist
The bp list format correspond to a bracket-list format readable in Mathematica. This format does not include any sequence information and, hence, is not suitable for sequence-dependent computations.

The format encodes a base-pair frame as:
```
{{origin_x,origin_y,origin_z}, {{d1x,d1y,d1z}, {d2x,d2y,d2z}, {d3x,d3y,d3z}}}
```
An input file in bp list format is composed by a list of base pair frames descibed in the x3DNAbp format above (each on a separate line).

Input file guidelines
- Always make sure the input file contains well-formatted numerical data (i.e., for example, no *^ symbols).
- When using x3DNA formats, make sure that every line which is not containing data starts with a #; in particular, x3DNA sometimes add information in the first lines that should be commented.
- When using a format encoding base-pair frames make sure your frames are orthogonalized.

## Command-line File Parsing: emDNA_parser
This tool can be used to convert data in one of the x3DNA format to any other format. In case data in a x3DNA format are converted to the bplist format, the sequence information are lost.

Input commands available: --x3DNA-bp-step-params-input=\<string\>, --x3DNA-bp-input=\<string\>, --bp-list-input=\<string\> <br>
Output commands: --get-x3DNA-params\>\<string\>, --get-x3DNA-bp\>\<string\>, --get-bp-list\>\<string\> <br>
Make sure that within each string you indicate desired file extensions, such as .par, .txt, or (for x3DNAbp) .dat<br>
  
Example: 
```
$ emDNA_parser --x3DNA-bp-input=test.dat --get-x3DNA-params>test.par
```
For additional assistance, use ```$ emDNA_parser --help```

## Command-line Forcefield Packaging: emDNA_ff_packager
emDNA can use an external force field (that is, a set of intrinsic step parameters and a set of force constant matrices). In order to use an external force field, the data needs to be packaged using the tool emDNA_force_field.

### Preparing the data
The intrinsic step parameters and force constant matrices need to be described in two different text files with a specific syntax before being packaged.

The intrinsic step parameters file has to be formatted as:
```
AA={0, 0, 34.2857, 0, 0, 3.4}
AC={0, 0, 34.2857, 0, 0, 3.4}
AG={0, 0, 34.2857, 0, 0, 3.4}
AT={0, 0, 34.2857, 0, 0, 3.4}
CA={0, 0, 34.2857, 0, 0, 3.4}
CC={0, 0, 34.2857, 0, 0, 3.4}
CG={0, 0, 34.2857, 0, 0, 3.4}
CT={0, 0, 34.2857, 0, 0, 3.4}
GA={0, 0, 34.2857, 0, 0, 3.4}
GC={0, 0, 34.2857, 0, 0, 3.4}
GG={0, 0, 34.2857, 0, 0, 3.4}
GT={0, 0, 34.2857, 0, 0, 3.4}
TA={0, 0, 34.2857, 0, 0, 3.4}
TC={0, 0, 34.2857, 0, 0, 3.4}
TG={0, 0, 34.2857, 0, 0, 3.4}
TT={0, 0, 34.2857, 0, 0, 3.4}
```
The file has to contain the 16 possible combinations of step sequences (the order is not relevant) and for each combination a bracket-list vector of dimension 6 is required.

The file describing the force constant matrices follows the same rules, except that, the bracket-list vector has to be of dimension 36(=6*6).

For these two files the order of the data in the bracket-list is: tilt, roll, twist, shift, slide, and rise. For the force constant matrices it means that the flatten bracket-list vector corresponds to:
```
{tilt/tilt, tilt/roll, tilt/twist, tilt/shift, tilt/slide, tilt/rise, roll/tilt, roll/roll, roll/twist, ...}
```
### Packaging the force field
To package the force field one has to use the emDNA_force_field. The tool requires the two files described above and a model name. The model name will be used to load the force field in emDNA. This tool creates two output files based on the ways the force field is packaged: a binary version to be used with emDNA (.ff) and a JSON-formatted text version for checking purposes (i.e., verify that the packaged force field is equivalent to the data contained in the two input files).

  
A command line example:
```
$ emDNA_ff_packager --intrinsic-steps-input=steps.txt --force-constants-input=fmat.txt --model-name=my_force_field
```

For additional assistance, use ```$ emDNA_ff_packager --help```


  
## Command-line Optimizations: emDNA

The emDNA cli tool is what conducts the optimization calculation. This tool requires:
- An input file (either in the x3DNAparams, x3DNAbp, or bp-list format)
- A forcefield (either pre-loaded or an external forcefield)
- Specific end conditions to be applied to the initial base pair collection

Additional options can be place in the command line including a customized output name (without this, the output files will be name "emDNA_minim"). 
emDNA-cli creates two output files. The first is a logfile (described below). The second contains the optimized base pair collection in a file format that matches the input file used and will always have "_opt" added to the output file name.

### Input files
An input file containing the initial base-pair collection has to be provided using one of the three following arguments depending on the file format:

* ```--x3DNA-bp-input=<filename>```: for input file with x3DNAbp format,
* ```--x3DNA-bp-step-params-input=<filename>```: for input file with x3DNAparams format,
* ```--bp-list-input=<filename>```: for input file with bplist format.

### Forcefields ("sequence-dependent models")
The sequence-dependence model for the minimization has to be indicated by either referring to an implemented model or by providing an external force field.

* ```--DNA-seqdep-model=<string>```: string should be one of the implemented force fields (see list of force fields),
* ```--DNA-external-model=<filename (binary)>```: specifies a file containing a packaged force field (see external force field).

The list of --DNA-seqdep-models currently set within emDNA include:
- **IdealDNA**: Not sequence-dependent, this force field corresponds to a straight B-DNA-like rest state with an helical repeat of 10.5 pb. The force constants are isotropic and quasi-inextensible (persistence lengths for bending is 47.7 nm and 66.6 nm for twisting).
- **IdealDNA_304**: Not sequence-dependent, same as IdealDNA but with an intrinsic twist of 34.3 degrees (the helical repeat is still close to 10.5 bp).
- **AnisoDNA**: Not sequence-dependent, this force field corresponds to a straight B-DNA-like rest state with an helical repeat of 10.5 pb. The force constants are anisotropic and quasi-inextensible. The persistence lengths are the same as for IdealDNA.
- **AnisoDNA_304**: Not sequence-dependent, same as AnisoDNA but with an intrinsic twist of 34.3 degrees (the helical repeat is still close to 10.5 bp).
- **Olson1998**: Fully sequence-dependent force field from values based on Olson, et al, PNAS 1998, vol 95, number 19, 11163-11168.
  
### End conditions
The end conditions applied to the base-pair collection has to be indicated by using one of the following flags:

* ```--free-collection```: the collection is free of any imposed constraint (used for debug purposes mainly),
* ```--hold-last-origin```: the origin of the last base pair is imposed and cannot change through the minimization process (imposed end-to-end vector),
* ```--hold-last-bp```: the origin and orientation of the last base pair is imposed and cannot change through the minimization process (imposed end-to-end vector and end-to-end rotation).

### Options
Below is the list of all options that can be used to customize the optimization:

* ```--output-progress=<integer>```: frequency for writing the current configuration (during the minimization) to tmp_confs.txt (default=0, which disables output)
* ```--energy-progress```: displays the current energy during the minimization
* ```--quiet```: suppresses console output (useful when using emDNA in scripts)
* ```--minim-settings=<{vector}>```: changes the minimization settings; the vector should be of the form {max-iterations,dx,f,g,max-step-size}
* ```--frozen-steps=<{list}>```: specify the ranges of base-pair steps to be considered as frozen. The list should be of the form: step_i:step_j,step_k:step_l, which specifies the ranges [i,j] and [k,l]. The numbering of steps does not start at 0 and the range of steps is an inclusive list
* ```--output-name=<string>```: prefix name for output files (default="emDNA_minim")


A command line example:
```
$ emDNA --x3DNA-bp-step-params-input=test.par --DNA-seqdep-model=IdealDNA --energy-progress --hold-last-bp --output-name=test
```

For additional assistance, use ```$ emDNA --help```

  
## Command-line Protein Binding: emDNA_probind
Details coming soon.

## Command-line Check Collisions: emDNA_check_collisions
Coming Soon
  
## Command-line Force Probe: emDNA_force_probe
Coming Soon
  
## Command-line Topology Data: emDNA_topology
A useful tool for both initial and optimized structures, higher-order information can be collected regarding the base pair collection's Linking Number, Twist, and Writhe. The only requirement is some input file.
Additional options include:
* ```--twist-density```: outputs the twist density along the collection of base pairs
* ```--virtual-last-bp```: a useful tool for optimizing circular structures, this assumes the last base pair is virtual
In addition, to specify an output name, use of '>' must follow the input file

A command-line example:
```
$ emDNA_topology --x3DNA-bp-input=test.dat>test_topology.txt --virtual-last-bp
```
For additional assistance, use ```$ emDNA_topology --help```
  
## Additional Details
### The Log File
The log file contains various information about the minimization process. By default the log file is named emDNA_minim.log (see --output-name option to change it). Here is an example:
```
log file timestamp: 2014-07-23 09:23:24


>> command line input
--x3DNA-bp-input=<filename>: 
--x3DNA-bp-step-params-input=<filename>: 
--bp-list-input=<filename>: minicircle_68bp.txt
--pulling-force=<{vector}>: {0,0,-1.2}
--DNA-seqdep-model=<string>: IdealDNA
--DNA-external-model=<filename (binary)>: 
--output-name=<string>: emDNA_minim
--frozen-steps=<{list}>: 
--minim-settings=<{vector}>: 
--output-progress=<integer>: 5
--free-collection: false
--hold-last-origin: false
--hold-last-bp: false
--pull-last-bp: true
--quiet: false
--energy-progress: true

>> initial bp collection
number of base pairs: 69
number of bp steps: 68
sequence: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
step parameters:
# SEQ      N/A        N/A        N/A        N/A        N/A        N/A        SHIFT      SLIDE      RISE       TILT       ROLL       TWIST
A-T        0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
A-T        0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    -0.00000   3.40000    2.78699    4.50115    31.76471

>> minimization results
minimization description: Alglib gradient-based | EEDR collection
initial energy: 53.5908021757
final energy: 6.3190724345
angular gradient norm: 0.0036884676
distance gradient norm: 0.0000649606
# iterations: 1480
return code: EPSG

>> energy contributions for step parameter pairs
{1.3684335288, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000}
{0.0000000000, 1.4338175002, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000}
{0.0000000000, 0.0000000000, 0.0000000034, 0.0000000000, 0.0000000000, 0.0000000000}
{0.0000000000, 0.0000000000, 0.0000000000, 0.0004527981, 0.0000000000, 0.0000000000}
{0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0004645536, 0.0000000000}
{0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0005388812}
```

The various sections correspond to:
- ```timestamp```: indicates the date at which the minimization was launch
- ```command line input```: list all command-line options with their parsed values (useful to diagnose an option problem)
- ```initial bp collection```: contains information about the input and lists the initial bair pair step parameters in the x3DNAbp format
- ```minimization results```: contains the results of the minimization (description, initial and final energies, gradient norms, ...)
- ```return code```:
    - ```EPSG```: the minimization successfully reduced the gradient norm below the threshold
    - ```EPSX```: the minimization stopped because the optimal point was not changing more than the threshold
    - ```EPSF```: the minimization stopped because the objective function (energy) was not decreasing more than the threshold
    - ```MAXIT```: the minimization reached the maximum number of iterations (can be relaunched with the output as new input)
    - ```BADCONDS```: the conditions are not sufficient to minimize the energy (this should not happen)
    - ```BADGRAD```: an error occurred in the computation of the gradient (this should not happen)
    - ```UNKNOWN```: an error occurred and most likely there is a bug in the software or a bad input
In case the return code is ```EPSX``` or ```EPSF```, a new minimization can be performed using the output as new input and changing the minimization settings (see ```--minim-settings```).

- ```energy contributions for step parameter pairs```: This 6x6 matrix gives the contributions to the energy of the different modes of deformation; column and row order is similiar: ```{'tilt','roll','twist','shift','slide','rise'}```


