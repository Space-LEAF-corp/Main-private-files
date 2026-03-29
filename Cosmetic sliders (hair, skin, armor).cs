public class AppearanceTemplate
{
    public float DefaultSkinTone;
    public float DefaultHairStyle;
    public float DefaultHairColor;
    public float DefaultArmorStyle;
}

public class AppearanceData
{
    public float SkinTone;     // 0–1
    public float HairStyle;    // index or blend
    public float HairColor;    // 0–1 or HSV param
    public float ArmorStyle;   // 0–1
}
