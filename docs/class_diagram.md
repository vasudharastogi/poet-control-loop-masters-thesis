# POET Class Diagram

```mermaid
classDiagram
    class RuntimeParameters {
        +bool print_progress
        +uint32_t work_package_size
        +bool use_dht
        +uint32_t dht_size
        +uint32_t dht_snaps
        +bool use_interp
        +uint32_t interp_size
        +uint32_t interp_min_entries
        +uint32_t interp_bucket_entries
        +bool use_ai_surrogate
        +bool as_rds
        +bool as_qs
        +string out_ext
        +string out_dir
        +vector~double~ timesteps
        +uint32_t checkpoint_interval
        +uint32_t stab_interval
        +double zero_abs
        +vector~double~ mape_threshold
        +vector~uint32_t~ ctrl_cell_ids
        +Rcpp::List init_params
    }

    class Field {
        +GetProps() vector~string~
        +AsVector() vector~double~
        +GetRequestedVecSize() size_t
        +update(Field) void
        +asSEXP() SEXP
        +operator[](string) vector~double~
    }

    class InitialList {
        -RInside& R
        +InitialList(RInside&)
        +importList(Rcpp::List, bool) void
        +getChemistryInit() ChemistryInit
        +getDiffusionInit() DiffusionInit
        +getInitialGrid() Field
    }

    class ChemistryModule {
        +ChemistryModule(uint32_t, ChemistryInit, MPI_Comm)
        +simulate(double) void
        +getField() Field&
        +WorkerLoop() void
        +masterSetField(Field) void
        +masterEnableSurrogates(SurrogateSetup) void
        +SetControlCellIds(vector~uint32_t~) void
        +SetControlModule(ControlModule*) void
        +setProgressBarPrintout(bool) void
        +set_ai_surrogate_validity_vector(SEXP) void
        +MasterLoopBreak() void
        +GetChemistryTime() double
        +GetMasterLoopTime() double
        +GetWorkerIdleTimings() vector~double~
        +GetWorkerPhreeqcTimings() vector~double~
        +GetWorkerDHTHits() vector~uint64_t~
        +GetWorkerDHTEvictions() vector~uint64_t~
        -Field field
        -uint32_t work_package_size
        -MPI_Comm comm
    }

    class DiffusionModule {
        +DiffusionModule(DiffusionInit, Field)
        +simulate(double) void
        +getField() Field&
        +getTransportTime() double
        -Field field
    }

    class RInsidePOET {
        +getInstance()$ RInsidePOET&
        +parseEval(string) SEXP
        +parseEvalQ(string) void
        +operator[](string) Proxy
    }

    class ChemistryInit {
        +dht_species SpeciesList
        +ai_surrogate_input_script string
    }

    class DiffusionInit {
    }

    class SurrogateSetup {
        +vector~string~ species_names
        +array~double,2~ base_totals
        +bool has_id
        +bool use_dht
        +uint32_t dht_size
        +uint32_t dht_snaps
        +string out_dir
        +bool use_interp
        +uint32_t interp_bucket_entries
        +uint32_t interp_size
        +uint32_t interp_min_entries
        +bool use_ai_surrogate
    }

    class Main {
        +main(int, char**) int
        -parseInitValues(int, char**, RuntimeParameters&) int
        -init_global_functions(RInside&) void
        -call_master_iter_end(RInside&, Field&, Field&) void
        -RunMasterLoop(RInsidePOET&, RuntimeParameters&, DiffusionModule&, ChemistryModule&, ControlModule&) Rcpp::List
        -getControlCellIds(vector~uint32_t~&, int, MPI_Comm) void
        -getSpeciesNames(Field&&, int, MPI_Comm) vector~string~
        -getBaseTotals(Field&&, int, MPI_Comm) array~double,2~
        -getHasID(Field&&, int, MPI_Comm) bool
    }

    Main --> RuntimeParameters : uses
    Main --> InitialList : creates
    Main --> ChemistryModule : creates
    Main --> DiffusionModule : creates
    Main --> RInsidePOET : uses
    Main --> Field : exchanges
    
    InitialList --> RInsidePOET : uses
    InitialList --> Field : creates
    InitialList --> ChemistryInit : provides
    InitialList --> DiffusionInit : provides
    
    ChemistryModule --> Field : manages
    ChemistryModule --> ChemistryInit : initialized with
    ChemistryModule --> SurrogateSetup : configured with
    
    DiffusionModule --> Field : manages
    DiffusionModule --> DiffusionInit : initialized with
    
    ChemistryModule ..> DiffusionModule : exchanges Field data
    DiffusionModule ..> ChemistryModule : exchanges Field data
    
    RuntimeParameters --> ChemistryInit : contains
```

## Key Relationships

- **Main** orchestrates the entire simulation, coordinating between modules
- **InitialList** parses R configuration and initializes all modules
- **ChemistryModule** and **DiffusionModule** exchange data via **Field** objects
- **Field** is the core data structure representing the simulation grid
- **RInsidePOET** provides the R runtime interface (singleton pattern)
- **RuntimeParameters** holds all command-line and configuration parameters
- **SurrogateSetup** configures advanced features (DHT, interpolation, AI surrogate)

## Module Communication Flow

1. Main reads configuration via `parseInitValues()`
2. `InitialList` imports R scripts and creates initial `Field`
3. `ChemistryModule` and `DiffusionModule` are initialized with their respective configurations
4. In simulation loop:
   - `DiffusionModule.simulate()` updates transport field
   - `ChemistryModule` receives updated field via `update()`
   - `ChemistryModule.simulate()` computes chemistry
   - `DiffusionModule` receives updated field back
5. MPI communication handled internally by modules
