# local_search

## Usage:
```
local_search [OPTION]...
```

## Description:
Performing local search algorithm to find the greatest deviation of provided trajectory frames by calculating RMSD value between them.

## Options:
`-c CONFIG`                           provide config parameters via CONFIG file
`-h, --help`                          display this message and exit

## Parameters if no config provided:
In descriptions: `[type:default]` format is used,
where type is the type of parameter, and default is its default value.

Parameter|Type:Default|Description
-|-|-
`--trajectory=TRAJECTORY`             | `[string:]` | `[mandatory]` trajectory filename in .pdb format
`--time-limit=TIME`                   | `[double:1.0]` | max time in minutes for whole local search to finish
`--omp-threads=NUM`                   | `[double:0]` | omp threads number per one cpu core
`--write-as-csv=[true/false]`         | `[bool:false]` | each run of a program generates one line in CSV format
`--repetitions=REPS`                  | `[int:2]` | number of program executions
`--jump-chance=PROB`                  | `[double:0.1]` | probability of jumping from local area
`--random-frame-chance=PROB`          | `[double:0.01]` | probability of choosing random frame while swapping allocations
`--memory-size=SIZE`                  | `[double:0.1]` | [0, 1] where 0 is no memory, and 1 is remembering whole matrix
`--random-seed=[true/false]`          | `[bool:true]` | random seed for srand()
`--matrix-size=SIZE`                  | `[int:-1]` | limiting matrix to SIZE by SIZE, if -1 then SIZE is max for current trajectory file
`--show-logs=[true/false]`            | `[bool:true]` | show any logs in the console
`--show-rmsd-counter=[true/false]`    | `[bool:false]` | show rsmd counter in the console
`--show-current-best=[true/false]`    | `[bool:true]` | show current best value, works only if --show-logs is set
`--show-route-best=[true/false]`      | `[bool:false]` | show current route best value, works only if --show-logs is set


## Examples:
```
local_search -c config.yml
```
```
local_search --trajectory=traj.pdb --time-limit=0.5 --repetitions=5
```

### All bool possible values:
- maps to true:  `true`  `t` `1` `yes` `y` `on`  ` ` &larr; ( nothing, e.g. `--write-as-csv` )
- maps to false: `false` `f` `0` `no`  `n` `off`
