from msrcsim.statistics import summarize_replicates


def test_prevalence_summary():
    rows = [
        {"accepted":True,"terminal_pattern":"1010","persisted_to_target":True,"is_2_2_pattern":True,"q2_minus_q3":0.3,"dominant_topology":1},
        {"accepted":True,"terminal_pattern":"0000","persisted_to_target":False,"is_2_2_pattern":False,"q2_minus_q3":0.0,"dominant_topology":0},
    ]
    s = summarize_replicates(rows, 0.1)
    assert s["accepted_replicates"] == 2
    assert s["asymmetric_quartet_distribution"]["proportion"] == 0.5
