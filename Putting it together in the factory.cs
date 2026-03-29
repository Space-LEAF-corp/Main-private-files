public static class CharacterFactory
{
    public static PlayerCharacter CreateCharacter(
        CharacterTemplate template,
        CharacterCreationData data,
        LineageData lineage,
        PlayerProfile profile)
    {
        var ch = new PlayerCharacter();

        // Base from template + overrides
        ch.Name         = data.Name ?? template.Name;
        ch.Race         = data.Race ?? template.Race;
        ch.Class        = data.Class ?? template.Class;
        ch.Strength     = data.Strength ?? template.Strength;
        ch.Agility      = data.Agility ?? template.Agility;
        ch.Intelligence = data.Intelligence ?? template.Intelligence;

        // Synergy
        var synergy = SynergyDatabase.Get(ch.Race, ch.Class);
        if (synergy != null) ApplySynergy(ch, synergy);

        // Appearance
        ch.Appearance = data.Appearance ?? new AppearanceData {
            SkinTone    = template.DefaultAppearance.DefaultSkinTone,
            HairStyle   = template.DefaultAppearance.DefaultHairStyle,
            HairColor   = template.DefaultAppearance.DefaultHairColor,
            ArmorStyle  = template.DefaultAppearance.DefaultArmorStyle
        };

        // Skill tree
        ch.Lineage   = lineage;
        ch.SkillTree = SkillTreeGenerator.Generate(template, ch.Race, ch.Class, lineage);

        // Backstory
        ch.Backstory = BackstoryGenerator.Generate(ch.Race, ch.Class, lineage);

        ch.IsCustom = data.HasCustomizations();
        return ch;
    }
}
