import pytest
from mite_extras.processing.data_classes import (
    Changelog,
    ChangelogEntry,
    Entry,
    EnyzmeDatabaseIds,
    Enzyme,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionDatabaseIds,
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
def reaction_databaseid():
    return ReactionDatabaseIds(ec="EC 1.2.3.4", mite="MITE0000000", rhea="32647")


@pytest.fixture
def reactionex():
    return ReactionEx(
        substrate="CCC",
        products=["CCCO"],
        forbidden_products=["CCC"],
        isBalanced=False,
        isIntermediate=True,
        description="Nonexistent reaction",
    )


@pytest.fixture
def reactionsmarts():
    return ReactionSmarts(
        reactionSMARTS="[#6]-[#6]-[#6]>>[#6]-[#6]-[#6]-[#8]",
        isIterative=False,
    )


@pytest.fixture
def reaction(reactionsmarts, reactionex, evidence, reaction_databaseid):
    return Reaction(
        tailoring=["Hydrolysis"],
        description="A nonexistent reaction",
        reactionSMARTS=reactionsmarts,
        reactions=[reactionex],
        evidence=[evidence],
        databaseIds=reaction_databaseid,
    )


@pytest.fixture
def enzyme_databaseids():
    return EnyzmeDatabaseIds(genpept="AAM70353.1", mibig="BGC0000581")


@pytest.fixture
def enyzmeaux(enzyme_databaseids):
    return EnzymeAux(
        name="AbcD",
        description="A nonexistent auxiliary enzyme.",
        databaseIds=enzyme_databaseids,
    )


@pytest.fixture
def enzyme(enyzmeaux, enzyme_databaseids):
    return Enzyme(
        name="BcdE",
        description="A nonexisting enzyme.",
        databaseIds=enzyme_databaseids,
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


def test_evidence_to_html_valid(evidence):
    html_dict = evidence.to_html()
    assert html_dict["evidenceCode"] == ["Heterologous expression"]


def test_reaction_databaseid_to_json_valid(reaction_databaseid):
    json_dict = reaction_databaseid.to_json()
    assert json_dict["mite"] == "MITE0000000"


def test_reaction_databaseid_to_html_valid(reaction_databaseid):
    html_dict = reaction_databaseid.to_html()
    assert html_dict["mite"][0] == "MITE0000000"


def test_reactionex_to_json_valid(reactionex):
    json_dict = reactionex.to_json()
    assert json_dict["substrate"] == "CCC"


def test_reactionex_to_html_valid(reactionex):
    html_dict = reactionex.to_html()
    assert html_dict["substrate"][0] == "CCC"


def test_reactionsmarts_to_json_valid(reactionsmarts):
    json_dict = reactionsmarts.to_json()
    assert json_dict["reactionSMARTS"] == "[#6]-[#6]-[#6]>>[#6]-[#6]-[#6]-[#8]"


def test_reactionsmarts_to_html_valid(reactionsmarts):
    html_dict = reactionsmarts.to_html()
    assert html_dict["reactionSMARTS"][0] == "[#6]-[#6]-[#6]>>[#6]-[#6]-[#6]-[#8]"


def test_reaction_to_json_valid(reaction):
    json_dict = reaction.to_json()
    assert json_dict["tailoring"] == ["Hydrolysis"]


def test_reaction_to_html_valid(reaction):
    html_dict = reaction.to_html()
    assert html_dict["tailoring"] == ["Hydrolysis"]


def test_enzymedatabaseids_to_json_valid(enzyme_databaseids):
    json_dict = enzyme_databaseids.to_json()
    assert json_dict["uniprot"] == "Q8KND5"


def test_enzymedatabaseids_to_html_valid(enzyme_databaseids):
    html_dict = enzyme_databaseids.to_html()
    assert html_dict["mibig"][0] == "BGC0000581"


def test_enyzmeaux_to_json_valid(enyzmeaux):
    json_dict = enyzmeaux.to_json()
    assert json_dict["name"] == "AbcD"


def test_enyzmeaux_to_html_valid(enyzmeaux):
    html_dict = enyzmeaux.to_html()
    assert html_dict["name"] == "AbcD"


def test_enzyme_to_json_valid(enzyme):
    json_dict = enzyme.to_json()
    assert json_dict["name"] == "BcdE"


def test_enzyme_to_html_valid(enzyme):
    html_dict = enzyme.to_html()
    assert html_dict["name"] == "BcdE"


def test_changelogentry_to_json_valid(changelogentry):
    json_dict = changelogentry.to_json()
    assert json_dict["contributors"] == ["AAAAAAAAAAAAAAAAAAAAAAAA"]


def test_changelogentry_to_html_valid(changelogentry):
    html_dict = changelogentry.to_html()
    assert html_dict["contributors"] == ["AAAAAAAAAAAAAAAAAAAAAAAA"]


def test_changelog_to_json_valid(changelog):
    json_dict = changelog.to_json()
    assert json_dict["version"] == "next"


def test_changelog_to_html_valid(changelog):
    html_dict = changelog.to_html()
    assert html_dict["version"] == "next"


def test_entry_to_json_valid(entry):
    json_dict = entry.to_json()
    assert json_dict["status"] == "pending"


def test_entry_to_html_valid(entry):
    html_dict = entry.to_html()
    assert html_dict["status"] == "pending"
