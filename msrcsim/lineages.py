from dataclasses import dataclass, field
@dataclass
class Lineage:
    lineage_id:int; descendants:frozenset[str]; population_id:str; arrangement:int; node_newick:str; birth_time:float=0.0
@dataclass(frozen=True)
class GenealogyEvent:
    locus_id:int; event_index:int; backward_time:float; event_type:str; population_before:str; population_after:str
    lineage_ids_before:str; lineage_id_after:str; descendants_before:str; arrangements_before:str; arrangements_after:str
    event_rate:float|None; total_rate:float|None; frequency_A1:float|None; m01:float|None; m10:float|None; lambda0:float|None; lambda1:float|None
@dataclass(frozen=True)
class CoalescenceRecord:
    locus_id:int; coalescence_index:int; time:float; population:str; child_lineage_1:int; child_lineage_2:int; parent_lineage:int; arrangement:int; descendant_set_1:str; descendant_set_2:str
@dataclass
class GenealogyResult:
    locus_id:int; newick:str; topology_index:int; coalescence_times:tuple[float,...]; events:tuple[GenealogyEvent,...]; coalescences:tuple[CoalescenceRecord,...]
