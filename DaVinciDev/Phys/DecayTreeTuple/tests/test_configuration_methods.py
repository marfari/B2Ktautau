###############################################################################
# (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
import random, string

from Configurables import DecayTreeTuple, TupleToolDecay
from DecayTreeTuple import Configuration


def _random_string(length=10):
    """Return a string of randomly chosen lowercase characters."""
    return ''.join(
        random.choice(string.ascii_lowercase) for i in range(length))


def _create_decaytreetuple(name=''):
    """Return a DecayTreeTuple instance for test cases.

    Because the DecayTreeTuple object is a Configurable, instances are
    'global', which could mess up a test case if the DecayTreeTuple instance
    was created previously. Avoid this by adding a random string to the name.
    """
    return DecayTreeTuple(name + _random_string())


def test_add_branches():
    """The addBranches method should add members on the DTT that are instances
    of TupleToolDecay, with a name equal to the branch name.
    """
    # The decay descriptors don't need to be valid for this test
    branches = dict(A='', B='', C='')

    dtt = _create_decaytreetuple()
    dtt.addBranches(branches)

    for bname in branches:
        try:
            bobj = getattr(dtt, bname)
        except AttributeError:
            assert False, 'Attribute for branch {0!r} not found'.format(bname)

        assert isinstance(bobj, TupleToolDecay)


def test_add_branches_no_whitespace_branch_name():
    """addBranches should not allow whitespace in the branch name."""
    for char in string.whitespace:
        branches = {char: ''}

        dtt = _create_decaytreetuple()
        try:
            dtt.addBranches(branches)
            assert False, 'Was able to make a branch with whitespace character {0!r}'.format(
                char)
        except NameError:
            pass


def test_set_descriptor_template():
    """The setDescriptorTemplate method should create the correct set of
    branches and the correct decay descriptor.
    """
    descriptor_template = '${Dst}[D*(2010)+ -> ${D0}(D0 -> ${D0_K}K- ${D0_pi}pi+) ${Dst_pi}pi+]CC'

    # The head does not need to marked, it always is implicitly
    expected_decay = '[D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+]CC'
    expected_branches = {
        'Dst': '[D*(2010)+ -> (D0 -> K- pi+) pi+]CC',
        'Dst_pi': '[D*(2010)+ -> (D0 -> K- pi+) ^pi+]CC',
        'D0': '[D*(2010)+ -> ^(D0 -> K- pi+) pi+]CC',
        'D0_K': '[D*(2010)+ -> (D0 -> ^K- pi+) pi+]CC',
        'D0_pi': '[D*(2010)+ -> (D0 -> K- ^pi+) pi+]CC'
    }

    dtt = _create_decaytreetuple()
    dtt.setDescriptorTemplate(descriptor_template)
    assert dtt.Decay == expected_decay
    for bname, bdecay in dtt.Branches.items():
        assert bname in expected_branches, '{0!r} branch expected, not found'.format(
            bname)
        exp_bdecay = expected_branches[bname]
        assert bdecay == exp_bdecay, '{0!r} expected, found {1!r}'.format(
            exp_bdecay, bdecay)


def test_set_descriptor_template_complex():
    """The setDescriptorTemplate method should create the correct set of
    branches and the correct decay descriptor, given a (more) complex template.
    """
    # Captures both the Cabibbo-favoured and doubly Cabibbo-suppressed decays
    descriptor_template = '[${Dst}(D*(2010)+ -> ${D0}(D0 -> ${D0_K}K- ${D0_pi}pi+) ${Dst_pi}pi+), ${Dst}(D*(2010)+ -> ${D0}(D0 -> ${D0_K}K+ ${D0_pi}pi-) ${Dst_pi}pi+)]CC'

    # The head does not need to marked, it always is implicitly
    expected_decay = '[(D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+), (D*(2010)+ -> ^(D0 -> ^K+ ^pi-) ^pi+)]CC'
    expected_branches = {
        'Dst':
        '[(D*(2010)+ -> (D0 -> K- pi+) pi+), (D*(2010)+ -> (D0 -> K+ pi-) pi+)]CC',
        'Dst_pi':
        '[(D*(2010)+ -> (D0 -> K- pi+) ^pi+), (D*(2010)+ -> (D0 -> K+ pi-) ^pi+)]CC',
        'D0':
        '[(D*(2010)+ -> ^(D0 -> K- pi+) pi+), (D*(2010)+ -> ^(D0 -> K+ pi-) pi+)]CC',
        'D0_K':
        '[(D*(2010)+ -> (D0 -> ^K- pi+) pi+), (D*(2010)+ -> (D0 -> ^K+ pi-) pi+)]CC',
        'D0_pi':
        '[(D*(2010)+ -> (D0 -> K- ^pi+) pi+), (D*(2010)+ -> (D0 -> K+ ^pi-) pi+)]CC'
    }

    dtt = _create_decaytreetuple()
    dtt.setDescriptorTemplate(descriptor_template)
    assert dtt.Decay == expected_decay
    for bname, bdecay in dtt.Branches.items():
        assert bname in expected_branches, '{0!r} branch expected, not found'.format(
            bname)
        exp_bdecay = expected_branches[bname]
        assert bdecay == exp_bdecay, '{0!r} expected, found {1!r}'.format(
            exp_bdecay, bdecay)


if __name__ == '__main__':
    test_add_branches()
    test_add_branches_no_whitespace_branch_name()

    test_set_descriptor_template()
    test_set_descriptor_template_complex()
