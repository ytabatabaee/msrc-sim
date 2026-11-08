class NeutralBackground:
    name = "neutral"

    def initial_state(self, population, segment):
        return "neutral"

    def coalescence_rate(self, population, state_a, state_b, effective_size, time):
        return 1.0 / (2.0 * effective_size)

    def recombination_states(self, population, state, left, right, time):
        return state, state
