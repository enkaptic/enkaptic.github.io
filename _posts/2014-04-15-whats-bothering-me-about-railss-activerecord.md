---
title: What's bothering me about Rails's ActiveRecord
tags:
- Computing
---
{% highlight ruby %}
class Product < ActiveRecord::Base
  # All of my attributes are public
end
{% endhighlight %}

I expect there are very good reasons to ditch the whole idea of <a href="http://en.wikipedia.org/wiki/Encapsulation_(object-oriented_programming)">encapsulation</a> for [Ruby on Rails models](http://guides.rubyonrails.org/active_record_basics.html). I'd love to know what they are. But, for now, it just doesn't feel right.
