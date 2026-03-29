public static SkillNode[] GenerateSkillTree(CharacterTemplate template, Race race, ClassType cls, LineageData lineage)
{
    var tree = template.SkillTreeTemplate
                       .Select(n => CloneNode(n))
                       .ToArray();

    // Example: unlock a hidden node if Jarvondis lineage + HH class
    if (lineage.IsJarvondis && cls == ClassType.HiddenHuman)
    {
        var node = tree.FirstOrDefault(n => n.Id == "HH_JARVONDIS_SIGIL");
        if (node != null) node.Unlocked = true;
    }

    return tree;
}
