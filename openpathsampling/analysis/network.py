import logging
import itertools
import random

import pandas as pd

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject
import openpathsampling.volume
import openpathsampling.ensemble

logger = logging.getLogger(__name__)

def index_to_string(index):
    n_underscore = index / 26
    letter_value = index % 26
    mystr = "_"*n_underscore + chr(65+letter_value)
    return mystr

class TransitionNetwork(StorableNamedObject):
    """
    Subclasses of TransitionNetwork are the main way to set up calculations

    Attributes
    ----------
    sampling_ensembles
    all_ensembles
    sampling_transitions
    """
    def __init__(self):
        super(TransitionNetwork, self).__init__()

    @property
    def sampling_ensembles(self):
        """
        Ensembles from the sampling transitions, excluding special ensembles.
        """
        return sum([t.ensembles for t in self.sampling_transitions], [])

    @property
    def analysis_ensembles(self):
        """
        Ensembles from the analysis transitions, excluding special ensembles.
        """
        return sum([t.ensembles for t in self.transitions.values()], [])

    @property
    def all_ensembles(self):
        """
        All ensembles in the sampling transitions, including special
        ensembles.
        """
        all_ens = self.sampling_ensembles
        for special_dict in self.special_ensembles.values():
            all_ens.extend(special_dict.keys())
        return all_ens

    @property
    def sampling_transitions(self):
        """The transitions used in sampling"""
        try:
            return self._sampling_transitions
        except AttributeError:
            return None


class GeneralizedTPSNetwork(TransitionNetwork):
    """General class for TPS-based method.

    The main differences between fixed-length and flexible-length TPS is a
    small change in the ensemble. In implementation, this means that they
    use different transition classes, and that they have slightly different
    function signatures (fixed-length requires a length argument).

    To simplify this, and to make the docstrings readable, we make each
    class into a simple subclass of this GeneralizedTPSNetwork, which acts
    as an abstract class that manages most of the relevant code.

    Parameters
    ----------
    initial_states : list of :class:`.Volume`
        acceptable initial states
    final_states : list of :class:`.Volume`
        acceptable final states
    allow_self_transitions : bool
        whether self-transitions (A->A) are allowed; default is False

    Attributes
    ----------
    TransitionType : :class:`paths.Transition`
        Type of transition used here. Sets, for example, fixed or flexible
        pathlengths.
    """
    TransitionType = NotImplemented
    def __init__(self, initial_states, final_states,
                 allow_self_transitions=False, **kwargs):
        # **kwargs gets passed to the transition
        super(GeneralizedTPSNetwork, self).__init__()
        try:
            iter(initial_states)
        except TypeError:
            initial_states = [initial_states]
        try:
            iter(final_states)
        except TypeError:
            final_states = [final_states]

        self.special_ensembles = {None : {}}

        self.initial_states = initial_states
        self.final_states = final_states

        all_initial = paths.join_volumes(initial_states)
        if len(initial_states) > 1:
            all_initial.name = "|".join([v.name for v in initial_states])

        if set(initial_states) == set(final_states) or len(final_states) == 1:
            all_final = all_initial
        else:
            all_final = paths.join_volumes(final_states)
            all_final.name = "|".join([v.name for v in final_states])
        
        self._sampling_transitions = []
        for my_initial in initial_states:
            my_final_states = [final for final in final_states
                               if my_initial != final or allow_self_transitions]
            my_final = paths.join_volumes(my_final_states)
            if len(my_final_states) > 1:
                my_final.name = "|".join([v.name for v in my_final_states])
            if  len(self._sampling_transitions) == 0:
                self._sampling_transitions = [
                    self.TransitionType(my_initial, my_final, **kwargs)
                ]
            elif len(self._sampling_transitions) == 1:
                self._sampling_transitions[0].add_transition(my_initial, 
                                                             my_final)
            else:
                raise RuntimeError("More than one sampling transition for TPS?")

        self.transitions = {
            (initial, final) : self.TransitionType(initial, final, **kwargs)
            for (initial, final) in itertools.product(initial_states,
                                                      final_states)
            if initial != final
        }


    def to_dict(self):
        ret_dict = {
            'transitions' : self.transitions,
            'x_sampling_transitions' : self._sampling_transitions,
        }
        return ret_dict

    @property
    def all_states(self):
        """list of all initial and final states"""
        return list(set(self.initial_states + self.final_states))

    @classmethod
    def from_dict(cls, dct):
        network = cls.__new__(cls)
        super(GeneralizedTPSNetwork, network).__init__()
        network._sampling_transitions = dct['x_sampling_transitions']
        network.transitions = dct['transitions']
        return network

    @classmethod
    def from_state_pairs(cls, state_pairs, **kwargs):
        sampling = []
        transitions = {}
        initial_states = []
        final_states = []
        for (initial, final) in state_pairs:
            initial_states += [initial]
            final_states += [final]
            if len(sampling) == 1:
                sampling[0].add_transition(initial, final)
            elif len(sampling) == 0:
                sampling = [cls.TransitionType(initial, final, **kwargs)]
            else:
                raise RuntimeError("More than one sampling transition for TPS?")

            transitions[(initial, final)] = cls.TransitionType(initial, final,
                                                               **kwargs)
        
        dict_result = {
            'x_sampling_transitions' : sampling,
            'transitions' : transitions
        }
        dict_result.update(kwargs)
        network = cls.from_dict(dict_result)
        network.initial_states = initial_states
        network.final_states = final_states
        return network


    @classmethod
    def from_states_all_to_all(cls, states, allow_self_transitions=False,
                               **kwargs):
        return cls(states, states,
                   allow_self_transitions=allow_self_transitions, **kwargs)


class TPSNetwork(GeneralizedTPSNetwork):
    """
    Class for flexible pathlength TPS networks (2-state or multiple state).
    """
    TransitionType = paths.TPSTransition
    # we implement these functions entirely to fix the signature (super's
    # version allow arbitrary kwargs) so the documentation can read them.
    def __init__(self, initial_states, final_states,
                 allow_self_transitions=False):
        super(TPSNetwork, self).__init__(initial_states, final_states,
                                         allow_self_transitions)

    @classmethod
    def from_state_pairs(cls, state_pairs, allow_self_transitions=False):
        return super(TPSNetwork, cls).from_state_pairs(state_pairs)

    @classmethod
    def from_states_all_to_all(cls, states, allow_self_transitions=False):
        return super(TPSNetwork, cls).from_states_all_to_all(
            states, allow_self_transitions
        )


class FixedLengthTPSNetwork(GeneralizedTPSNetwork):
    """
    Class for fixed pathlength TPS networks (2-states or multiple states).
    """
    TransitionType = paths.FixedLengthTPSTransition
    # as with TPSNetwork, we don't really need to add these functions.
    # However, without them, we need to explicitly name `length` as
    # length=value in these functions. This frees us of that, and gives us
    # clearer documentation.
    def __init__(self, initial_states, final_states, length,
                 allow_self_transitions=False):
        super(FixedLengthTPSNetwork, self).__init__(
            initial_states, final_states,
            allow_self_transitions=allow_self_transitions, length=length
        )

    @classmethod
    def from_state_pairs(cls, state_pairs, length):
        return super(FixedLengthTPSNetwork, cls).from_state_pairs(
            state_pairs, length=length
        )

    @classmethod
    def from_states_all_to_all(cls, states, length,
                               allow_self_transitions=False):
        return super(FixedLengthTPSNetwork, cls).from_states_all_to_all(
            states=states,
            allow_self_transitions=allow_self_transitions,
            length=length
        )


class TISNetwork(TransitionNetwork):
    # NOTE: this is an abstract class with several properties used by many
    # TIS-based networks
    # TODO: most of the analysis stuff should end up in here; the bigger
    # differences are in setup, not analysis
    def __init__(self):
        super(TISNetwork, self).__init__()

    def from_transitions(self, transitions, interfaces=None):
        # this will have to be disabled until I can do something
        # better with it
        pass

    def set_fluxes(self, flux_dictionary):
        """
        Parameters
        ----------
        flux_dictionary : dict of 2-tuple to float
            keys are in the form (state, interface), and values are the
            associated flux

        Raises
        ------
        KeyError
            If the flux for one of the transitions isn't in the dictionary.
        """
        # for now, if you don't have all the fluxes needed, it raises a
        # KeyError
        for trans in self.transitions.values():
            trans._flux = flux_dictionary[(trans.stateA, trans.interfaces[0])]


    @property
    def minus_ensembles(self):
        return self.special_ensembles['minus'].keys()

    @property
    def ms_outers(self):
        return self.special_ensembles['ms_outer'].keys()

    @property
    def all_states(self):
        return list(set(self.initial_states + self.final_states))

    def get_state(self, snapshot):
        """
        Find which core state a snapshot is in, if any

        Parameters
        ----------
        snapshot : `openpathsampling.engines.BaseSnapshot`
            the snapshot to be tested

        Returns
        -------
        `openpathsampling.Volume`
            the volume object defining the state
        """
        for state in self.all_states:
            if state(snapshot):
                return state

        return None


#def join_mis_minus(minuses):
    #pass

#def msouter_state_switching(mstis, steps):

class MSTISNetwork(TISNetwork):
    """
    Multiple state transition interface sampling network.

    The way this works is that it sees two effective sets of transitions.
    First, there are sampling transitions. These are based on ensembles
    which go to any final state. Second, there are analysis transitions.
    These are based on ensembles which go to a specific final state.

    Sampling is done using the sampling transitions. Sampling transitions
    are stored in the `from_state[state]` dictionary. For MSTIS, the flux
    and total crossing probabilities are independent of the final state, and
    so the analysis calculates them in the sampling transitions, and copies
    the results into the analysis transitions. This way flux and total
    crossing probably are only calculated once per interface set.

    The conditional transition probability depends on the final state, so it
    (and the rate) are calculated using the analysis transitions. The
    analysis transitions are obtained using `.transition[(stateA, stateB)]`.
    """
    def to_dict(self):
        ret_dict = { 
            'from_state' : self.from_state,
            'states' : self.states,
            'special_ensembles' : self.special_ensembles,
            'trans_info' : self.trans_info
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        network = cls.__new__(cls)

        # replace automatically created attributes with stored ones
        network.from_state = dct['from_state']
        network.special_ensembles = dct['special_ensembles']
        network.states = dct['states']
        network.__init__(
            trans_info=dct['trans_info']
        )
        return network

    def __init__(self, trans_info):
        """
        Creates MSTISNetwork, including interfaces.

        Parameters
        ----------
        trans_info : list of tuple
            Details of each state-based ensemble set. 3-tuple in the order
            (state, interfaces, orderparameter) where state is a Volume,
            interfaces is a list of Volumes, and orderparameters is a
            CollectiveVariable
        """
        super(MSTISNetwork, self).__init__()
        self.trans_info = trans_info
        # build sampling transitions
        if not hasattr(self, "from_state"):
            self.special_ensembles = {}
            self.from_state = {}
            self.build_fromstate_transitions(trans_info)

        self._sampling_transitions = self.from_state.values()

        # by default, we set assign these values to all ensembles
        self.hist_args = {}

        self.transitions = {}
        self.build_analysis_transitions()

    @property
    def all_states(self):
        return self.states

    def build_analysis_transitions(self):
        # set up analysis transitions (not to be saved)
        for stateA in self.from_state.keys():
            state_index = self.states.index(stateA)
            fromA = self.from_state[stateA]
            other_states = self.states[:state_index]+self.states[state_index+1:]
            for stateB in other_states:
                trans = paths.TISTransition(
                    stateA=stateA,
                    stateB=stateB,
                    interfaces=fromA.interfaces,
                    name=str(stateA) + "->" + str(stateB),
                    orderparameter=fromA.orderparameter
                )
                # override created stuff
                trans.ensembles = fromA.ensembles
                trans.minus_ensemble = fromA.minus_ensemble
                self.transitions[(stateA, stateB)] = trans

#    def disallow(self, stateA, stateB):

    def build_fromstate_transitions(self, trans_info):
        """
        Builds the sampling transitions (the self.from_state dictionary).

        This also sets self.states (list of states volumes), self.outers
        (list of interface volumes making the MS-outer interface), and 
        self.outer_ensembles (list of TISEnsembles associated with the
        self.outers interfaces). Additionally, it gives default names
        volumes, interfaces, and transitions.

        Parameters
        ----------
        trans_info : list of 4-tuples
            See description in __init__.

        """
        states, interfaces, orderparams = zip(*trans_info)
        # NAMING STATES (give default names)
        all_states = paths.volume.join_volumes(states).named("all states")
        all_names = list(set([s.name for s in states]))
        unnamed_states = [s for s in states if not s.is_named]
        name_index = 0
        for state in unnamed_states:
            while index_to_string(name_index) in all_names:
                name_index += 1
            state.named(index_to_string(name_index))
            name_index += 1

        # BUILDING ENSEMBLES
        outer_ensembles = []
        self.states = states
        for (state, ifaces, op) in trans_info:
            state_index = states.index(state)
            other_states = states[:state_index]+states[state_index+1:]
            union_others = paths.volume.join_volumes(other_states)
            union_others.named("all states except " + str(state.name))

            this_trans = paths.TISTransition(
                stateA=state, 
                stateB=union_others,
                interfaces=ifaces[:-1],
                name="Out " + state.name,
                orderparameter=op
            )

            self.from_state[state] = this_trans

            this_minus = self.from_state[state].minus_ensemble
            this_inner = self.from_state[state].ensembles[0]
            try:
                self.special_ensembles['minus'][this_minus] = [this_trans]
            except KeyError:
                self.special_ensembles['minus'] = {this_minus : [this_trans]}


            outer_ensemble = paths.TISEnsemble(
                initial_states=state,
                final_states=all_states,
                interface=ifaces[-1]
            )
            outer_ensemble.named("outer " + str(state))
            outer_ensembles.append(outer_ensemble)

        ms_outer = paths.ensemble.join_ensembles(outer_ensembles)
        transition_outers = self.from_state.values()
        try:
            self.special_ensembles['ms_outer'][ms_outer] = transition_outers
        except KeyError:
            self.special_ensembles['ms_outer'] = {ms_outer : transition_outers}


    def __str__(self):
        mystr = "Multiple State TIS Network:\n"
        for state in self.from_state.keys():
            mystr += str(self.from_state[state])
        return mystr


    def rate_matrix(self, steps, force=False):
        """
        Calculate the matrix of all rates.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            steps to be analyzed
        force : bool (False)
            if True, cached results are overwritten

        Returns
        -------
        pandas.DataFrame
            Rates from row_label to column_label. Diagonal is NaN.
        """
        # for each transition in from_state:
        # 1. Calculate the flux and the TCP
        self._rate_matrix = pd.DataFrame(columns=self.states,
                                         index=self.states)
        for stateA in self.from_state.keys():
            transition = self.from_state[stateA]
            # set up the hist_args if necessary
            for histname in self.hist_args.keys():
                trans_hist = transition.ensemble_histogram_info[histname]
                if trans_hist.hist_args == {}:
                    trans_hist.hist_args = self.hist_args[histname]
        
            transition.total_crossing_probability(steps=steps,
                                                  force=force)
            transition.minus_move_flux(steps=steps, force=force)
            for stateB in self.from_state.keys():
                if stateA != stateB:
                    analysis_trans = self.transitions[(stateA, stateB)]
                    analysis_trans.copy_analysis_from(transition)


        for trans in self.transitions.values():
            rate = trans.rate(steps)
            self._rate_matrix.set_value(trans.stateA, trans.stateB, rate)
            #print trans.stateA.name, trans.stateB.name, 
            #print rate

        return self._rate_matrix

    def generate_inter_state_samples(self, snapshot, engine):
        vol_all = reduce(lambda x, y: x | y, self.states)
        states_missing = list(self.states)  # we want a copy

        # all inter-core trajectories
        core2core = paths.SequentialEnsemble(
            [
                paths.OptionalEnsemble(paths.AllInXEnsemble(vol_all)),
                paths.AllOutXEnsemble(vol_all),
                paths.SingleFrameEnsemble(
                    paths.AllInXEnsemble(vol_all)
                )
            ]
        )

        # core X -> core X trajectories
        core2core_list = {
            state: paths.SequentialEnsemble([
                paths.SingleFrameEnsemble(
                    paths.AllInXEnsemble(state)),
                paths.AllOutXEnsemble(state),
                paths.SingleFrameEnsemble(
                    paths.AllInXEnsemble(state))
            ]) for state in self.states
        }

        # start at template
        last_traj = []
        count = 0

        start_snap = snapshot
        samp_list = []
        start_snap_state = {state: [] for state in self.states}
        total_length = 0

        paths.tools.refresh_output(
            '       TOTAL FRAMES [%6d] // missing: [ %s ]\n'
            '----------------------------\n'
            '%s' % (
                total_length,
                ' '.join([s.name for s in states_missing]),
                '\n'.join(last_traj[-10:])
            )
        )


        all_states_found = False

        while not all_states_found:
            new_traj_part = engine.generate_forward(
                start_snap,
                core2core
            )

            # look through all missing states
            part_ensemble = None
            to_delete = []
            for idx, state in enumerate(states_missing):
                if core2core_list[state](new_traj_part):
                    # found X -> X path so stop looking for state X
                    part_ensemble = core2core_list[state]
                    to_delete = [state]
                    break

            states_missing = sorted(list(set(states_missing) - set(to_delete)))

            last_traj.append(
                '[%4d] %3s -> %3s   [%6d]' % (
                    count,
                    self.get_state(new_traj_part[0]).name
                    if self.get_state(new_traj_part[0]) else '',
                    self.get_state(new_traj_part[-1]).name,
                    len(new_traj_part)
                )
            )

            samp_list.append(
                paths.Sample(
                    replica=0,
                    trajectory=new_traj_part,
                    ensemble=part_ensemble,
                )
            )

            all_states_found = not states_missing

            # the last frame as a potential starting point
            last_frame = new_traj_part[-1].reversed
            start_snap_state[self.get_state(last_frame)].append(last_frame)

            # pick a starting point either from the last position or
            # from a state that has been visited but where we have no X->X path
            starting_points = \
                [last_frame] + \
                sum([points for state, points in start_snap_state.items()
                     if state in states_missing], [])

            start_snap = random.choice(starting_points)

            total_length += len(new_traj_part)

            # update output
            paths.tools.refresh_output(
                '       TOTAL FRAMES [%6d] // missing: [ %s ]\n'
                '----------------------------\n'
                '%s' % (
                    total_length,
                    ' '.join([s.name for s in states_missing]),
                    '\n'.join(last_traj[-10:])
                )
            )

            count += 1

        return samp_list

    def generate_initial_sampleset(self, samples, engine):
        return paths.SampleSet.generate_from_sampleset(
            self.all_ensembles,
            samples,
            engine
        )

#def multiple_set_minus_switching(mistis, steps):

class MISTISNetwork(TISNetwork):
    """
    Multiple interface set TIS network.

    Input is given as a list of 4-tuples. Each 4-tuple represents a
    transition, and is in the order: 
        (initial_state, interfaces, order_parameter, final_states)
    This will create the `input_transitions` objects.

    Attributes
    ----------
    input_transitions : list of TISTransition
        the transitions given as input
    sampling_transitions : list of TISTransition
        the transitions used in sampling
    transitions : list of TISTransition
        the transitions used in analysis

    Note
    ----
        The distinction between the three types of transitions in the object
        are a bit subtle, but important. The `input_transitions` are, of
        course, the transitions given in the input. These are A->B
        transitions, but would allow any other state. The
        `sampling_transitions` are what are used in sampling. These are
        A->any transitions if strict sampling is off, or "A->B & not_others"
        if strict sampling is on. Finally, the regular `transitions` are the
        transitions that are used for analysis (use the sampling ensembles
        for the interfaces, but also A->B).

    Parameters
    ----------
    trans_info : list of tuple
        Details of each interface set. 4-tuple in the order (initial_state,
        interfaces, orderparameter, final_state) where initial_state and
        final_state are Volumes, interfaces is a list of Volumes, and
        orderparameter is a CollectiveVariable
    strict_sampling : bool
        whether the final state from the tuple is the *only* allowed final
        state in the sampling; default False
    """
    # NOTE: input_transitions are in addition to the sampling_transitions
    # and the transitions (analysis transitions)
    def __init__(self, trans_info, strict_sampling=False):
        super(MISTISNetwork, self).__init__()
        self.trans_info = trans_info
        self.strict_sampling = strict_sampling
        states_A, interfaces, orderparams, states_B = zip(*trans_info)
        self.initial_states = list(set(states_A))
        self.final_states = list(set(states_B))
        list_all_states = list(set(self.initial_states + self.final_states))

        # name states
        all_state_names = list(set([s.name for s in list_all_states]))
        unnamed_states = [s for s in list_all_states if not s.is_named]
        name_index = 0
        for state in unnamed_states:
            while index_to_string(name_index) in all_state_names:
                name_index += 1
            state.named(index_to_string(name_index))
            name_index += 1


        if not hasattr(self, "input_transitions"):
            self.input_transitions = {
                (stateA, stateB) :
                paths.TISTransition(stateA, stateB, interface, orderparam,
                                    name=stateA.name+"->"+stateB.name)
                for (stateA, interface, orderparam, stateB) in self.trans_info
            }

        if not hasattr(self, 'x_sampling_transitions'):
            self.special_ensembles = {}
            self.build_sampling_transitions(self.input_transitions.values())
        self._sampling_transitions = self.x_sampling_transitions


        # by default, we set assign these values to all ensembles
        self.hist_args = {}

        self.build_analysis_transitions()


    def to_dict(self):
        ret_dict = {
            'special_ensembles' : self.special_ensembles,
            'transition_pairs' : self.transition_pairs,
            'x_sampling_transitions' : self.x_sampling_transitions,
            'transition_to_sampling' : self.transition_to_sampling,
            'input_transitions' : self.input_transitions,
            'trans_info' : self.trans_info,
            'strict_sampling' : self.strict_sampling
        }
        return ret_dict

    @staticmethod
    def from_dict(dct):
        network = MISTISNetwork.__new__(MISTISNetwork)
        network.special_ensembles = dct['special_ensembles']
        network.transition_pairs = dct['transition_pairs']
        network.transition_to_sampling = dct['transition_to_sampling']
        network.input_transitions = dct['input_transitions']
        network.x_sampling_transitions = dct['x_sampling_transitions']
        network.__init__(dct['trans_info'], dct['strict_sampling'])
        return network


    def build_sampling_transitions(self, transitions):
        # identify transition pairs
        for initial in self.initial_states:
            transition_pair_dict = {}
            for t1 in [t for t in transitions if t.stateA==initial]:
                reverse_trans = None
                for t2 in transitions:
                    if t2.stateA==t1.stateB and t2.stateB==t1.stateA:
                        transition_pair_dict[t1] = t2
            # TODO: speed this up with a set?
            for key in transition_pair_dict.keys():
                value = transition_pair_dict[key]
                if value in transition_pair_dict.keys():
                    del transition_pair_dict[value]
        self.transition_pairs = [(k, transition_pair_dict[k]) 
                                 for k in transition_pair_dict.keys()]

        all_in_pairs = reduce(list.__add__, map(lambda x: list(x), 
                                                self.transition_pairs))

        # build sampling transitions
        all_states = paths.join_volumes(self.initial_states + self.final_states)
        all_states_set = set(self.initial_states + self.final_states)
        self.transition_to_sampling = {}
        for transition in transitions:
            stateA = transition.stateA
            stateB = transition.stateB
            if self.strict_sampling:
                final_state = stateB
                other_states = paths.join_volumes(all_states_set -
                                                  set([stateA, stateB]))
                ensemble_to_intersect = paths.AllOutXEnsemble(other_states)
            else:
                final_state = all_states
                ensemble_to_intersect = paths.FullEnsemble()
            # TODO: fix following for strict_sampling
            if transition not in all_in_pairs:
                # if we don't have a pair partner, use all interfaces
                sample_trans = paths.TISTransition(
                    stateA=stateA,
                    stateB=final_state,
                    interfaces=transition.interfaces,
                    orderparameter=transition.orderparameter
                )
            else:
                # if we do have a pair partner, outermost is MS-interface
                sample_trans = paths.TISTransition(
                    stateA=stateA,
                    stateB=final_state,
                    interfaces=transition.interfaces[:-1],
                    orderparameter=transition.orderparameter
                )
            new_ensembles = [e & ensemble_to_intersect 
                             for e in sample_trans.ensembles]
            if self.strict_sampling:
                for (old, new) in zip(new_ensembles, sample_trans.ensembles):
                    old.name = new.name + " strict"
            sample_trans.ensembles = new_ensembles
            sample_trans.named("Sampling " + str(stateA) + "->" + str(stateB))
            self.transition_to_sampling[transition] = sample_trans

        self.x_sampling_transitions = self.transition_to_sampling.values()

        # build non-transition interfaces 

        # combining the MS-outer interfaces
        for pair in self.transition_pairs:
            this_outer = paths.ensemble.join_ensembles(
                [pair[0].ensembles[-1], pair[1].ensembles[-1]]
            )
            s_pair = [self.transition_to_sampling[p] for p in pair]
            try:
                self.special_ensembles['ms_outer'][this_outer] = list(s_pair)
            except KeyError:
                self.special_ensembles['ms_outer'] = {this_outer : list(s_pair)}

        
        # combining the minus interfaces
        for initial in self.initial_states:
            innermosts = []
            trans_from_initial = [
                t for t in self.x_sampling_transitions
                if t.stateA==initial
            ]
            for t1 in trans_from_initial:
                innermosts.append(t1.interfaces[0])
            minus = paths.MinusInterfaceEnsemble(
                state_vol=initial,
                innermost_vols=innermosts
            )
            try:
                self.special_ensembles['minus'][minus] = trans_from_initial
            except KeyError:
                self.special_ensembles['minus'] = {minus : trans_from_initial}

    def build_analysis_transitions(self):
        self.transitions = {}
        for trans in self.input_transitions.values():
            sample_trans = self.transition_to_sampling[trans]
            stateA = trans.stateA
            stateB = trans.stateB
            analysis_trans = paths.TISTransition(
                stateA=stateA,
                stateB=stateB,
                interfaces=sample_trans.interfaces,
                orderparameter=sample_trans.orderparameter
            )
            analysis_trans.ensembles = sample_trans.ensembles
            analysis_trans.named(trans.name)
            #analysis_trans.special_ensembles = sample_trans.special_ensembles
            self.transitions[(stateA, stateB)] = analysis_trans


    def rate_matrix(self, steps, force=False):
        self._rate_matrix = pd.DataFrame(columns=self.final_states,
                                         index=self.initial_states)
        for trans in self.transitions.values():
            # set up the hist_args if necessary
            for histname in self.hist_args.keys():
                trans_hist = trans.ensemble_histogram_info[histname]
                if trans_hist.hist_args == {}:
                    trans_hist.hist_args = self.hist_args[histname]
            tcp = trans.total_crossing_probability(steps=steps,
                                                   force=force)
            if trans._flux is None:
                logger.warning("No flux for transition " + str(trans.name)
                               + ": Rate will be NaN")
                trans._flux = float("nan")
                # we give NaN so we can calculate the condition transition
                # probability automatically

            rate = trans.rate(steps)
            self._rate_matrix.set_value(trans.stateA, trans.stateB, rate)

        return self._rate_matrix
