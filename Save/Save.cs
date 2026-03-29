public static string SerializeCharacter(PlayerCharacter ch)
{
    var save = MapToSave(ch);
    return JsonUtility.ToJson(save); // or your JSON lib
}
