public class CrossGameShadowBank
{
    private readonly ShadowInventory _shadow;
    private readonly CrossGameEconomy _economy;

    public CrossGameShadowBank(ShadowInventory shadow, CrossGameEconomy economy)
    {
        _shadow = shadow;
        _economy = economy;
    }

    public bool ConvertAndStore(GlobalItemId fromId, GlobalItemId toId, long quantity, IItem toItemDef)
    {
        string fromKey = fromId.ToString();
        string toKey = toId.ToString();

        if (_shadow.GetQuantity(fromKey) < quantity)
            return false;

        if (!_economy.TryConvert(fromKey, toKey, quantity, out long converted))
            return false;

        _shadow.RemoveItem(fromKey, (int)quantity);

        // Here, toItemDef would be a definition bound to toKey
        return _shadow.AddItem(toItemDef, (int)converted);
    }
}
