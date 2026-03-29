public static class CharacterFactory
{
    public static PlayerCharacter CreateCharacter(CharacterTemplate template, CharacterCreationData data)
    {
        PlayerCharacter character = new PlayerCharacter();

        // Start with template defaults
        character.Name = data.Name ?? template.Name;
        character.Class = data.Class ?? template.Class;
        character.Strength = data.Strength ?? template.Strength;
        character.Agility = data.Agility ?? template.Agility;
        character.Intelligence = data.Intelligence ?? template.Intelligence;

        character.IsCustom = data.HasCustomizations();

        return character;
    }
}
