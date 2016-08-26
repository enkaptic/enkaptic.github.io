---
title: WordPress to Rails (5) static pages
tags:
- Computing
redirect_to: http://www.anthonysmith.me.uk/2013/10/21/wordpress-to-rails-5-static-pages/
---
WordPress lumps posts and static pages together in one database table, but they need to be handled quite differently.

At least for the time being, I'm leaving them in the same table.

First, in the <code>Post</code> model, a different permalink is needed for static pages and for posts:

{% highlight ruby %}
def permalink
  if post_type == 'page'
    post_name
  else
    post_date.strftime("%Y/%m/%d") + (post_name.empty? ? '' : "/#{post_name}")
  end
end
{% endhighlight %}

Then, also in the <code>Post</code> model, I have several scopes:

{% highlight ruby %}
scope :published_posts, where(:post_status => 'publish', :post_type => 'post').order('post_date DESC')
scope :published_pages, where(:post_status => 'publish', :post_type => 'page')
scope :pages, where(:post_type => 'page')
scope :not_pages, where('post_type != page')
{% endhighlight %}

Then two different methods to find posts, depending on whether the user is logged in or not (see below):

{% highlight ruby %}
def self.find_by_permalink(permalink)
  if permalink !~ /\//
    pages.find(:first, :conditions => ["post_name = ?", permalink])
  else
    year, month, day, name = permalink.split('/')
    name ||= ''
    not_pages.find(:first, :conditions => ["YEAR(post_date) = ? AND MONTH(post_date) = ? AND DAYOFMONTH(post_date) = ? AND post_name = ?", year.to_i, month.to_i, day.to_i, name])
  end
end

def self.find_published_by_permalink(permalink)
  if permalink !~ /\//
    published_pages.find(:first, :conditions => ["post_name = ?", permalink])
  else
    year, month, day, name = permalink.split('/')
    name ||= ''
    published_posts.find(:first, :conditions => ["YEAR(post_date) = ? AND MONTH(post_date) = ? AND DAYOFMONTH(post_date) = ? AND post_name = ?", year.to_i, month.to_i, day.to_i, name])
  end
end
{% endhighlight %}

Then, in <code>PostController</code>, the <code>show</code> method becomes a bit more advanced (or messy):

{% highlight ruby %}
def show
  # Show all posts if logged in, otherwise restrict to published posts
  if User.find_by_ID(session[:user_id])
    @post = Post.find_by_permalink(params[:id]) || not_found
  else
    @post = Post.find_published_by_permalink(params[:id])
    if not @post
      # Check if post _would_ be visible if logged in; ask to log in if so
      if Post.find_by_permalink(params[:id])
        session[:return_to] = request.url
        redirect_to login_url, :notice => "Please log in"
        return
      else
        not_found
      end
    end
  end
  ...
end
{% endhighlight %}

The <code>not_found</code> method is in <code>ApplicationController</code> (see <a href="http://stackoverflow.com/questions/2385799/how-to-redirect-to-a-404-in-rails/7099193">here</a>):

{% highlight ruby %}
def not_found
  raise ActionController::RoutingError.new('Not Found')
end
{% endhighlight %}

<code>routes.rb</code> needs to be changed, to route (say) <code>.../about</code> to the appropriate page, with the following line going <em>after</em> all the other routing information (but before the <code>root</code> statement):

{% highlight ruby %}
resources :posts, :only => [:show, :update, :destroy, :edit], :path => '', :id => /\d{4}\/\d{2}\/\d{2}(\/[0-9a-z\-_]+)?|[0-9a-z\-_]+/
{% endhighlight %}

And there are some changes to the post view, to hide details that are not relevant to static pages (such as date and comments).
