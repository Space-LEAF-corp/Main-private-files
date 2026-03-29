public static PlayerCharacter DeserializeCharacter(string json, CharacterTemplate baseTemplate)
{
    var save = JsonUtility.FromJson<PlayerCharacterSave>(json);
    var ch = new PlayerCharacter();
    // map back + optionally re-apply template defaults if fields missing
    return ch;
}
