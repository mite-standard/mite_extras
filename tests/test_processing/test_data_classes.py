import pytest
from mite_extras.processing.data_classes import (
    Changelog,
    ChangelogEntry,
    Entry,
    Enzyme,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionEx,
    ReactionSmarts,
)


@pytest.fixture
def evidence():
    return Evidence(
        evidenceCode=["Heterologous expression"],
        references=["doi:10.1002/cbic.201200016"],
    )


@pytest.fixture
def reactionex(evidence):
    return ReactionEx(
        substrate="CCC",
        products=["CCCO", "OCCC"],
        isBalanced=False,
        isIntermediate=True,
        description="Nonexistent reaction",
        databaseIds=["MITE0000000"],
        evidence=[evidence],
    )


@pytest.fixture
def reactionsmarts(evidence):
    return ReactionSmarts(
        reactionSMARTS="[#6]-[#6]>>[#6]-[#6]-[#8]",
        isIterative=False,
        databaseIds=["MITE0000000"],
        evidence=[evidence],
    )


@pytest.fixture
def reaction(reactionsmarts, reactionex):
    return Reaction(
        tailoring=["Hydrolysis"],
        description="A nonexistent reaction",
        reactionSMARTS=reactionsmarts,
        reactions=[reactionex],
    )


@pytest.fixture
def enyzmeaux():
    return EnzymeAux(
        name="AbcD",
        description="A nonexistent auxiliary enzyme.",
        databaseIds=["uniprot:Q9X2V9"],
    )


@pytest.fixture
def enzyme(enyzmeaux):
    return Enzyme(
        name="BcdE",
        description="A nonexisting enzyme.",
        databaseIds=["uniprot:Q9X2V8"],
        auxiliaryEnzymes=[enyzmeaux],
        references=["doi:10.1128/jb.181.8.2659-2662.1999"],
    )


@pytest.fixture
def changelogentry():
    return ChangelogEntry(
        contributors=["AAAAAAAAAAAAAAAAAAAAAAAA"],
        reviewers=["AAAAAAAAAAAAAAAAAAAAAAAA"],
        date="0000-00-00",
        comment="An example comment",
    )


@pytest.fixture
def changelog(changelogentry):
    return Changelog(version="next", date="0000-00-00", entries=[changelogentry])


@pytest.fixture
def entry(changelog, enzyme, reaction):
    return Entry(
        accession="MITE0000000",
        quality="high",
        status="pending",
        retirementReasons=["Example reason"],
        changelog=[changelog],
        enzyme=enzyme,
        reactions=[reaction],
        comment="An example comment",
        attachments={"attachment1": "An example attachment"},
    )


def test_evidence_to_json_valid(evidence):
    json_dict = evidence.to_json()
    assert json_dict["evidenceCode"] == ["Heterologous expression"]


def test_reactionex_to_json_valid(reactionex):
    json_dict = reactionex.to_json()
    assert json_dict["substrate"] == "CCC"


def test_reactionsmarts_to_json_valid(reactionsmarts):
    json_dict = reactionsmarts.to_json()
    assert json_dict["reactionSMARTS"] == "[#6]-[#6]>>[#6]-[#6]-[#8]"


def test_reaction_to_json_valid(reaction):
    json_dict = reaction.to_json()
    assert json_dict["tailoring"] == ["Hydrolysis"]


def test_enyzmeaux_to_json_valid(enyzmeaux):
    json_dict = enyzmeaux.to_json()
    assert json_dict["name"] == "AbcD"


def test_enzyme_to_json_valid(enzyme):
    json_dict = enzyme.to_json()
    assert json_dict["name"] == "BcdE"


def test_changelogentry_to_json_valid(changelogentry):
    json_dict = changelogentry.to_json()
    assert json_dict["contributors"] == ["AAAAAAAAAAAAAAAAAAAAAAAA"]


def test_changelog_to_json_valid(changelog):
    json_dict = changelog.to_json()
    assert json_dict["version"] == "next"


def test_entry_to_json_valid(entry):
    json_dict = entry.to_json()
    assert json_dict["status"] == "pending"
